# 使用官方的 Ubuntu 基础镜像
FROM nvidia/cuda:11.2.2-devel-ubuntu20.04 AS builder-image

# 避免在安装时出现提示
ARG DEBIAN_FRONTEND=noninteractive

# 更新软件包列表
RUN apt-get update

# 安装 Java、curl、以及必要的工具
RUN apt-get install -y openjdk-11-jre-headless curl gnupg2 lsb-release software-properties-common

# 添加 Docker 的官方 GPG 密钥
RUN curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add -

# 设置 Docker 稳定版仓库
RUN add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"

# 安装 Docker CE（社区版）
RUN apt-get update && apt-get install -y docker-ce docker-ce-cli containerd.io

# 清理缓存
RUN apt-get clean && rm -rf /var/lib/apt/lists/*

# 安装 Nextflow
RUN curl -s https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/
    
RUN chmod -R 777 /mnt/

