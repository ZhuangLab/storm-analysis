"""Run the SCons build as part of the Python build.

The C libraries are built by SCons rather than by setuptools, so without
this they have to be built by hand before 'pip install .'. Getting that
order wrong is easy and the result is confusing: the install succeeds, and
then fails later when it tries to load a library that was never compiled.
Re-running SCons afterwards does not fix an already installed copy either.

See [tool.setuptools.cmdclass] in pyproject.toml for how this is wired in.
"""
import os
import platform
import subprocess
import sys


def runSCons():
    """
    Build the C libraries.

    The 'scons' script is not reliably on PATH inside pip's isolated build
    environment, so call the entry point it would have called. SCons itself
    comes from [build-system] requires.
    """
    args = [sys.executable,
            "-c", "from SCons.Script.Main import main; main()",
            "-Q"]

    # Windows has no default compiler the way the other platforms do, so
    # SConstruct has to be told. Set SA_SCONS_COMPILER to use something
    # other than mingw.
    compiler = os.environ.get("SA_SCONS_COMPILER", "")
    if (len(compiler) == 0) and (platform.system() == "Windows"):
        compiler = "mingw"
    if (len(compiler) > 0):
        args.append("compiler=" + compiler)

    print("Building C libraries:", " ".join(args), flush = True)
    subprocess.check_call(args, cwd = os.path.dirname(os.path.abspath(__file__)))


from setuptools.command.build_py import build_py as _build_py


class build_py(_build_py):

    def run(self):
        runSCons()
        super().run()
