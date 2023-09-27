# Exercise-22-02-2022
Exercise 3 on the course of Simulation of Biological Systems

The recording of the hands-on session can be found on the MoBioChem youtube channel, in the following link:

https://www.youtube.com/watch?v=NRpnKUJXA7E&ab_channel=MoBioChem

Please either git clone this repository before starting the exercise, or copy the 3_REDOX folder with all of its contents on your working dierctory.

# XTB Installation
Apart from the orca software, we will also use xtb. Therefore, please:
1) Download xtb from

https://github.com/grimme-lab/xtb/releases/tag/v6.4.1

depending on you operative system. It comes precompiled so you really just need to unzip the corresponding folder. Then unzip it:

```
tar -xf xtb-6.4.1-linux-x86_64.tar.xz
```

2) Add these variables to your .bashrc file (or alternatively, paste the following on your terminal):

```
export XTBHOME=/path/to/xtb/bin
export XTBPATH=/path/to/xtb/bin
export PATH=$XTBHOME:$PATH
export PATH=$XTBPATH:$PATH
```

where /path/to/xtb is the location where you unzipped xtb. For example, if you downloaded the version 6.4.1 and you put it on your home, it should be:

```
/home/pepito/xtb-6.4.1/bin
```

if your username is pepito.

3) Very Important: we need to put the xtb executable on the ORCA installation folder with the name **otool_xtb**. We can create a link on that folder as follows:

```
ln -sf /path/to/xtb/bin/xtb /path/to/orca/otool_xtb
```

where /path/to/xtb/bin is the same as before (notice the "second" xtb at the end: that is the actual xtb executable, which you are linking to on the orca folder).


# Regarding the Second Part of the Exercise
For the second part (2_qmmm_sp) we will use trajectories having 10 snapshots (instead of 100), so please consider that the range is 10 in that part of the exercise.
