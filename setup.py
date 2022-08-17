from setuptools import setup

version = '0.1.0'
install_requires=["Bio","numpy"]

setup(
    name='codonOptimizer',
    version=version,
    description="Tools to identify preferred codons and codon optimize coding sequences",
    keywords='codon',
    author='Alex Crocker',
    author_email='alex.w.crocker@gmail.com',
    python_requires='>=3.3',
    install_requires=install_requires,
)
