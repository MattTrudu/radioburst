from setuptools import setup

setup(
    name='radioburst',
    version='0.1.0',
    description='A Python package to manipulate Fast Radio Burst Data',
    url='https://github.com/MattTrudu/radioburst.git',
    author='Matt Trudu',
    author_email='trudumatteo@outlook.com',
    license='BSD 2-clause',
    packages=['radioburst'],
    install_requires=['your',
                      'numpy',
                      ],

    classifiers=[
        'Development Status :: 1 - Planning',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Operating System :: POSIX :: Linux',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
    ],
)
