from setuptools import setup

setup(name='scito',
      version='0.2.0',
      description='YeLab software to perform SCITO-seq analysis',
      url='https://github.com/yelabucsf/SCITO-seq',
      author='Anton Gvaihir Ogorodnikov, Ye Lab UCSF',
      author_email='anton.ogorodnikov@ucsf.edu',
      license='GNU V3',
      packages=['scito',],
      install_requires=[
            'numpy',
            'seaborn',
            'matplotlib',
            'pandas',
            'scikit-learn',
            'scipy',
            'statsmodels',
      ],
      scripts=[],
      zip_safe=False)
