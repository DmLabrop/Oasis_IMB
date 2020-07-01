#!/usr/bin/env python

from setuptools import setup

# Version number
major = 2020
minor = 1

setup(name = "oasis_imb",
      version = "%d.%d" % (major, minor),
      description = "Oasis_IMB - Immersed Boundary Navier-Stokes solvers in FEniCS based on the oasis solver \
      originally developed by developed by Mikael Mortensen et al..",
      classifiers = [
          'Intended Audience :: Developers',
          'Intended Audience :: Science/Research',
          'Programming Language :: Python ',
          'License :: OSI Approved :: GNU Library or Lesser General Public License (LGPL)',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Software Development :: Libraries :: Python Modules',
          ],
      packages = ["oasis",
                  "oasis.problems",
                  "oasis.problems.NSfracStep",
                  "oasis.problems.NSCoupled",
                  "oasis.solvers",
                  "oasis.solvers.NSfracStep",
                  "oasis.solvers.NSfracStep.LES",
                  "oasis.solvers.NSCoupled",
                  "oasis.common",
                  ],
      package_dir = {"oasis": "oasis"},
      entry_points = {'console_scripts': ['oasis=oasis.run_oasis:main']},
    )