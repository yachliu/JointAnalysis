from setuptools import setup, find_packages

setup(
    name='jointanalysis',
    version='0.1',
    packages=find_packages(),
    description='A short description of your project',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='louie',
    author_email='ychlouie@gmail.com',
    url='https://github.com/yachliu/JointAnalysis',
    install_requires=[
        # 依赖列表
        'sqlite == 3.36.0',
        'numba == 0.53.1',
        'PyMSNumpress == 0.2.3',
        'zlib == 1.2.11',
        'numpy == 1.19.2',
        'lmdb == 4.6.3',
        'pandas == 1.1.5',
        'xgboost == 1.3.3',
        'seaborn == 0.11.2',
        'matplotlib == 3.3.4',
        'pulearn == 0.0.7',
        'tqdm == 4.62.3',
        'loguru == 0.7.2',
        'click == 8.0.3',
        'scikit-learn == 0.24.2',
        'scipy == 1.5.2',
        'statsmodels == 0.12.2'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
)