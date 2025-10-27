from setuptools import setup


setup(
    name="genespace", 
    version='0.1.0',
    author="Konstantin Zaytsev",
    url='https://github.com/conzaytsev/GeneSpace',
    description='Gene Space',
    long_description='Python module for analysis of genomic sequences as points a multidimensional Gene Space.',
    packages=['genespace'],
    install_requires=['scikit-learn', 'numpy'],
    python_requires='>3.9',
    classifiers= [
        'License :: OSI Approved :: GPL-3.0 License',
        'Programming Language :: Python :: 3'
    ]
)
