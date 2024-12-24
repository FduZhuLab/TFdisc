from setuptools import setup, find_packages

setup(
    name='TFdisc',                  
<<<<<<< HEAD
    version='1.2.0',               
    description='in silico transcription factor (TF) perturbation simulator', 
    author='Haiyang Wang',             
    author_email='your.email@example.com',  
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
=======
    version='1.1.0',               
    description='in silico transcription factor (TF) perturbation simulator', 
    author='Haiyang Wang',             
    author_email='your.email@example.com',  
>>>>>>> 144b3617b0f9ac7af893c5570c0093b8ddc5b51a
    url='https://github.com/ocean-debug/TFdisc',
    packages=find_packages(),       
    install_requires=[              
        'numpy', 'pandas', 'scikit-learn', 'arboreto', 'tqdm', 'dask', 'rpy2','scanpy'
    ],
    classifiers=[                 
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',       
)

