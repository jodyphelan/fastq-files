import setuptools

version = [l.strip() for l in open("fastq_files/__init__.py") if "version" in l][0].split('"')[1]

setuptools.setup(

	name="fastq_files",

	version=version,
	packages=["fastq_files"],
	license="GPLv3",
	long_description="Toolkit to find fastq files",
	scripts= [
        'scripts/find-fastq-files.py'
		],
)
