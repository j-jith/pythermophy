from setuptools import setup

def readme():
    with open('README.rst') as fileh:
        return fileh.read()

setup(name = 'pythermophy',
      version = '0.1',
      description = 'Predict density, speed of sound and heat capacities of pure fluids',
      long_description=readme(),
      classifiers = [
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2.7',
        'Topic :: Scientific/Engineering'
      ],
      url = 'http://github.com/j-jith/pythermophy',
      author = 'Jithin Jith',
      author_email = 'j.jith@outlook.com',
      license = 'MIT',
      packages = ['pythermophy'],
      install_requires = [
          'numpy',
          'scipy'
      ],
      include_package_data=True,
      zip_safe = False)
