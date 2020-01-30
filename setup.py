import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='SOAPswga',
    version='0.0.1',
    author='Jane Yu',
    author_email='janeyu@berkeley.edu',
    classifiers=[
        'Programming Language :: Python',
        'Development Status :: 4 - Beta',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License (GPL)'
    ],  
    packages=setuptools.find_packages(),
    url='https://github.com/janeaisley/soapswga',
    license='LICENSE.txt',
    description='Pipeline to select primer sets for selective whole-genome amplification.',
    long_description=open('README.md').read(),
    entry_points={'console_scripts': ['soapswga = src.pipeline:main']}
)
