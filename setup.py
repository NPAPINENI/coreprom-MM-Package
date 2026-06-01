from setuptools import setup, find_packages

setup(
    name="coreprom-MM",
    version="0.1.0",
    packages=find_packages(),
    python_requires=">=3.8",
    install_requires=["numpy", "pandas", "scikit-learn", "rpy2"],
    entry_points={
        "console_scripts": [
            "coreprom-props=coreprom_MM.cli:main",
        ],
    },
)
