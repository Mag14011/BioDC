from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="BioDC",
    version="3.0.0",
    author="Matthew J. Guberman-Pfeffer and Caleb L. Herron",
    author_email="Matthew_Guberman-Pfe@baylor.edu",
    description="A Python program that automates and accelerates the computation of redox potentials, cooperativities, and conductivities in (polymeric) multi-heme cytochromes.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Mag14011/BioDC",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",  
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    include_package_data=True,  
    package_data={
        'biodc': ['data/forcefield/*.lib', 'data/forcefield/*.frcmod'],
    },
    entry_points={
        'console_scripts': [
            'biodc-cli=biodc.biodc_cli:main',            # comprehensive interactive commandline 
            'biodc-spr-gui=biodc.biodc_struc_prep_gui:main', # structure preparaiton gui 
            'biodc-ee-gui=biodc.biodc_eng_eval_gui:main',   # energetic evalulation gui (under development)
        ],
    },
)
