import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="orbipyd",
    version="1",
    author="Denis Zagorodnev",
    author_email="denis.zagorodnev@gmail.com",
    description="A package which must be my masterpice!",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/BoberSA/OrbiPy",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)