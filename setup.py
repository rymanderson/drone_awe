import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="drone_awe",
    version="1.0.1",
    author="Tyler Critchfield and Ryan Anderson",
    author_email="rymanderson@gmail.com",
    description="A versatile modeling package for drone performance in various weather environments.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    download_url="https://github.com/rymanderson/drone_awe.git",
    url="https://github.com/rymanderson/drone_awe",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    install_requires=[
        'scipy',
        'numpy',
        'gekko',
        'matplotlib',
        'datetime',
    ],
)
