from skbuild import setup

URL = "https://github.com/kjelljorner/libconeangle"
DESCRIPTION = "Python library for libconeangle"
LONG_DESCRIPTION = f"""\
{DESCRIPTION}. For more information, see the [project repository]({URL}).
"""

setup(
    name="libconeangle",
    version="0.1.1",
    author="Kjell Jorner",
    author_email="kjell.jorner@gmail.com",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    long_description_content_type="text/markdown",
    url=URL,
    packages=["libconeangle"],
    python_requires=">=3.8",
    install_requires=["numpy>=1.20"],
    include_package_data=True,
    cmake_args=["-DSKBUILD=ON"],
    license="MIT",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
)
