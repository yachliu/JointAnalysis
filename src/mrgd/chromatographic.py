import numba
import PyMSNumpress
import zlib
import numpy as np

from bisect import bisect

# def organize_chroms_dimension(run_p_chroms: list, run_p_rts: list) -> np.array:
#     pr_chrom = np.array(run_p_chroms[0])
#     fr_chroms = np.array(run_p_chroms[1:])
#     pr_rt = run_p_rts[0]
#     fr_rt = run_p_rts[1]

#     for left_i in range(len(fr_rt)):
#         left_rtid = bisect(pr_rt, fr_rt[left_i]) - 1
#         if left_rtid >= 0:
#             break
#     if (left_i - left_rtid) == 0:
#         pass
#     elif left_rtid - left_i > 0:
#         pr_rt = pr_rt[left_rtid - left_i :]
#         pr_chrom = pr_chrom[left_rtid - left_i :]
#     elif left_i - left_rtid > 0:
#         fr_chroms = fr_chroms[:, left_i - left_rtid:]
#         fr_rt = fr_rt[left_i - left_rtid:]

#     rev_fr_rt = list(-np.array(fr_rt)[::-1])
#     rev_pr_rt = list(-np.array(pr_rt)[::-1])

#     for right_i in range(len(rev_pr_rt)):
#         right_rtid = bisect(rev_fr_rt, rev_pr_rt[right_i]) - 1
#         if right_rtid >= 0:
#             break
#     if (right_i - right_rtid) == 0:
#         pass
#     elif right_rtid - right_i > 0:
#         fr_rt = fr_rt[: -(right_rtid - right_i)]
#         fr_chroms = fr_chroms[:, : -(right_rtid - right_i)]
#     elif right_i - right_rtid > 0:
#         pr_rt = pr_rt[: -(right_i - right_rtid)]
#         pr_chrom = pr_chrom[: -(right_i - right_rtid)]

#     run_p_xics = np.concatenate((pr_chrom[np.newaxis, :], fr_chroms), axis = 0)

#     return run_p_xics, pr_rt


def organize_chroms_dimension(run_p_chroms: list, run_p_rts: list) -> np.array:
    pr_chrom = np.array(run_p_chroms[0])
    fr_chroms = np.array(run_p_chroms[1:])
    pr_rt = np.array(run_p_rts[0])
    fr_rt = np.array(run_p_rts[1])

    if len(pr_rt) == len(fr_rt):
        return np.array(run_p_chroms), pr_rt
    else:

        for psi in range(len(pr_rt)):
            psi_left = pr_rt[psi]
            psi_right = pr_rt[psi+1]
            sleft_label = fr_rt >= psi_left
            sright_label = fr_rt < psi_right
            s_label = sleft_label & sright_label
            if np.sum(s_label) > 0:
                break
        fsi = np.where(s_label)[0][0]

        for fei in range(len(fr_rt) - 1, 0, -1):
            fei_left = fr_rt[fei - 1]
            fei_right = fr_rt[fei]
            sleft_label = pr_rt > fei_left
            sright_label = pr_rt <= fei_right
            s_label = sleft_label & sright_label
            if np.sum(s_label) > 0:
                break
        pei = np.where(s_label)[0][0]

        part_pr_chrom = pr_chrom[psi: pei+1]
        part_fr_chroms = fr_chroms[:, fsi:fei+1]
        part_pr_length = len(part_pr_chrom)
        part_fr_length = len(part_fr_chroms[0])
        run_pr_rt = pr_rt[psi: pei+1]

        if part_pr_length > part_fr_length:
            part_pr_chrom = part_pr_chrom[: part_fr_length]
            run_pr_rt = run_pr_rt[: part_fr_length]
        elif part_pr_length < part_fr_length:
            part_fr_chroms = part_fr_chroms[:, : part_pr_length]
        else:
            pass
        run_p_xics = np.concatenate((part_pr_chrom[np.newaxis, :], part_fr_chroms), axis = 0)

        assert run_p_xics.shape[1] == len(run_pr_rt), print("Error")

        return run_p_xics, run_pr_rt

@numba.jit(nopython = True)
def smooth_array_nb(arr):
    # arr = np.array(arr)
    new_arr = np.zeros_like(arr)
    new_arr[0] = (2 * arr[0] + arr[1]) / 3
    new_arr[-1] = (2 * arr[-1] + arr[-2]) / 3
    for x in range(1, len(arr) - 1):
        new_arr[x] = (0.5*arr[x] + 0.25*arr[x-1] + 0.25*arr[x+1])
    return new_arr

def get_chromid_from_rid_native_id(db, search_key):
    value = db.get(search_key.encode('utf-8'))
    return value.decode('utf-8') if value else None

def organize_chroms(p_u_rids, rid2chrom_conn, p_native_ids, rid_native2chromid_db, logger):
    
    p_rid2xics, p_rid2rts = {}, {}
    for p_u_rid in p_u_rids:
        run_conn = rid2chrom_conn[p_u_rid]
        run_cur = run_conn.cursor()
        run_p_rts = []
        run_p_chroms = []
        for p_native_id in p_native_ids:
            run_p_chrom_id = int(get_chromid_from_rid_native_id(rid_native2chromid_db[p_u_rid], str(p_native_id)))
            comp_data = run_cur.execute(f'SELECT COMPRESSION, DATA FROM DATA WHERE CHROMATOGRAM_ID == {run_p_chrom_id}').fetchall()

            for comp, xic in comp_data:
                result = []
                if comp == 5:
                    PyMSNumpress.decodeLinear(zlib.decompress(xic), result)
                    run_p_rts.append(result)
                elif comp == 6:
                    PyMSNumpress.decodeSlof(zlib.decompress(xic), result)
                    run_p_chroms.append(smooth_array_nb(np.array(result)).tolist())
                else:
                    logger.error("Error: PyMSNumpress")
                    raise
        run_cur.close()
        run_p_xics, run_pr_rt = organize_chroms_dimension(run_p_chroms, run_p_rts)
        p_rid2xics[p_u_rid] = run_p_xics
        p_rid2rts[p_u_rid] = run_pr_rt
    return p_rid2xics, p_rid2rts

def return_xics_by_fid(p_table, p_rid2xics, p_rid2rts, logger):
# 获取fid对应的xics
    p_fid2xics = {}
    p_length_xics = []
    for idx in p_table.index:
        fid, rid, s_rt, e_rt = p_table.loc[idx, ["ID", "RUN_ID", "LEFT_WIDTH", "RIGHT_WIDTH"]]
        t_rts = p_rid2rts[rid]
        t_xics = p_rid2xics[rid]

        s_rt_ind = bisect(t_rts, s_rt) - 1
        s_rt_ind = 0 if s_rt_ind < 0 else s_rt_ind
        e_rt_ind = bisect(t_rts, e_rt)
        e_rt_ind = len(t_rts) if e_rt_ind > len(t_rts) else e_rt_ind

        tt_xics = t_xics[:, s_rt_ind : e_rt_ind]
        p_fid2xics[fid] = tt_xics
        p_length_xics.append(tt_xics.shape[1])
    p_xic_unified_length = np.ceil(np.percentile(p_length_xics, 95)).astype(int)
    
    p_fid2xics_pad = {}
    for k, v in p_fid2xics.items():
        n_pad = p_xic_unified_length - v.shape[1]
        left_n_pad = n_pad // 2
        right_n_pad = n_pad - left_n_pad
        if left_n_pad == 0 and right_n_pad == 0:
            p_fid2xics_pad[k] = v
        elif left_n_pad >= 0 and right_n_pad >= 0:
            p_fid2xics_pad[k] = np.pad(v, ((0, 0), (right_n_pad, left_n_pad)), 'constant', constant_values=(0, 0))
        elif left_n_pad <= 0 and right_n_pad <= 0:
            p_fid2xics_pad[k] = v[:, -right_n_pad : left_n_pad]
        else:
            logger.error(f'Error: {right_n_pad}, {left_n_pad}')
    return p_fid2xics_pad, p_xic_unified_length