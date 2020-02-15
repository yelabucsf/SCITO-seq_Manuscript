from setuptools import setup

setup(name='scito',
      version='0.0.1',
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
            'anndata',
            'scanpy',
            'scikit-learn',
            'pyclustering',
      ],
      scripts=[],
      zip_safe=False)
