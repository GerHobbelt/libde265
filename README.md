
libde265 internals
==========================================

![libde265](libde265.png)

The internal branch offers a lot of more insight into the decoder, allowing to extract detailed CU/TU data and much more.
We use this version in YUView (https://github.com/IENT/YUView) to show statistics.

Building
========

Status of the internal branch:
[![Build Status](https://ci.appveyor.com/api/projects/status/github/ChristianFeldmann/libde265?branch=internal)](https://ci.appveyor.com/project/ChristianFeldmann/libde265)
[![Build Status](https://travis-ci.org/ChristianFeldmann/libde265.svg?branch=internal)](https://travis-ci.org/ChristianFeldmann/libde265)

The generated files can be found at https://dl.bintray.com/christianfeldmann/libde265_Internal/

If you got libde265 from the git repository, you will first need to run
the included `autogen.sh` script to generate the `configure` script.

libde265 has no dependencies on other libraries, but both optional example programs
have dependencies on:

- SDL (optional for dec265's YUV overlay output),

- Qt (required for sherlock265),

- libswscale (required for sherlock265 if libvideogfx is not available).

- libvideogfx (required for sherlock265 if libswscale is not available,
  optional for dec265).

Libvideogfx can be obtained from
  http://www.dirk-farin.net/software/libvideogfx/index.html
or
  http://github.com/farindk/libvideogfx


You can disable building of the example programs by running `./configure` with
<pre>
  --disable-dec265        Do not build the dec265 decoder program.
  --disable-sherlock265   Do not build the sherlock265 visual inspection program.
</pre>

Additional logging information can be turned on and off using these `./configure` flags:
<pre>
  --enable-log-error      turn on logging at error level (default=yes)
  --enable-log-info       turn on logging at info level (default=no)
  --enable-log-trace      turn on logging at trace level (default=no)
</pre>


Build using cmake
=================

cmake scripts to build libde265 and the sample scripts `dec265` and `enc265` are
included and can be compiled using these commands:

```
mkdir build
cd build
cmake ..
make
```

See the [cmake documentation](http://www.cmake.org) for further information on
using cmake on other platforms.


Prebuilt binaries
=================

Binary packages can be obtained from this [launchpad site](https://launchpad.net/~strukturag/+archive/libde265).


Software using libde265
=======================

Libde265 has been integrated into these applications:

- gstreamer plugin, [source](https://github.com/strukturag/gstreamer-libde265), [binary packages](https://launchpad.net/~strukturag/+archive/libde265).

- VLC plugin [source](https://github.com/strukturag/vlc-libde265), [binary packages](https://launchpad.net/~strukturag/+archive/libde265).

- Windows DirectShow filters, https://github.com/strukturag/LAVFilters/releases

- ffmpeg fork, https://github.com/farindk/ffmpeg

- ffmpeg decoder [source](https://github.com/strukturag/libde265-ffmpeg)

- libde265.js JavaScript decoder [source](https://github.com/strukturag/libde265.js), [demo](https://strukturag.github.io/libde265.js/).


License
=======

The library `libde265` is distributed under the terms of the GNU Lesser
General Public License. The sample applications are distributed under
the terms of the GNU General Public License.

See `COPYING` for more details.

Copyright (c) 2013-2014 Struktur AG
Contact: Dirk Farin <farin@struktur.de>
