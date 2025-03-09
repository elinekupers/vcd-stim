# vcd-stim

MATLAB code repository to create and present the stimuli of the Visual Cognition Dataset (VCD).


## Goal

Create and run the 39 different types of visual experiments of the Visual Cognition Dataset across 30 MRI sessions. Visual experiments are based on crossings of 5 stimulus classes and 9 task classes.

Stimulus classes:
* Gabors
* RDKs: Random dot motion kinetograms, 
* Dot: single simple dot, 
* Cobj: complex objects
* NS: natural scenes.

Task classes:
* FIX: fixation
* CD: contrast-detection
* SCC: superclass categorization
* PC: perceptual categorization
* WM: working memory
* LTM: long-term memory
* IMG: imagery
* WHAT: object categorization
* WHERE: object localization
* HOW: scene/object affordance).


## Dependencies

This stimulus presentation code is written for the iMac that lives permanently at the 7TAS MRI at the Center for Magnetic Resonance Research (CMRR) at the University of Minnesota. 

This code requires [Psychtoolbox-3](https://github.com/Psychtoolbox-3/Psychtoolbox-3).

The 7TAS iMac specs:
* Retina 5K, 27-inch, 2017
* macOS High Sierra (version 10.13.6), 
* Processor 4.2 GHz Intel Core i7
* Memory 16 GB 2400 MHz DDR4
* Graphics Radeon Pro 575 4096 MB
* MATLAB version 2016b and R2017b
* Psychtoolbox version 3.0.14 (December 30, 2016).

This code has also been tested on a MacbookPro with macOS 10.14.6 Mojave, and MATLAB version 2018b, PTB version (November 17, 2020), git commit ef093cbf296115badddb995fa06452e34c8c7d02.

The stimulus presentation code does not run on MacBook Pro's with Silicon chips (M1 or M2).

## Code and folder overview

* _runme_vcdcore.m_ 	This is the main function to run the core experiment of VCD.
* _vcd_rootPath.m_ 		Function to set the rootpath to relative to the base of this folder.

* external 				Folder with external functions from other toolboxes	
* functions				Folder with vcd related functions 
* ptb 					Folder with vcd functions that run stimulus presentation with psychtoolbox
* scripts				Folder with standalone scripts to create stimuli
* tinkers				Folder with code tinkering around and such (probably should be removed at some point)
* utils 				Folder with small and simple utility functions.


Ignored-by-git folders:
* workspaces			Folder where the stimuli and stimulus info files live.
* figs					Folder where debug figures are stored
* data 					Folder where subject's button presses and created stimuli are stored (if requested).


## Terminology

* _Monitor refresh rate:_ The rate by which a monitor runs vertical retrace (VBL). For example, BOLD screen monitor refresh is 120 Hz (for reference: 6 frames = (1000/120) * 6 = 50 ms). Most LCD monitors have a refresh rate of 60 Hz
* _Frames/ framerate / framedur:_ Because we present our stimuli at a slower rate than the monitor refresh rate, we “bundle” or “skip” several monitor refreshes into a “frame”. We run our stimuli at 30 Hz, so we bundle 4 120 Hz (BOLDscreen) monitor refreshes or 2 60 Hz (regular LCD) monitor refreshes.
* _Image:_  The retinal image: the image on the BOLD screen that falls onto the retina.
* _Session:_ The 2-3 hr scan session at the MRI scanner
* _Run:_ Each session consists of approx. 10 runs (+ 2 resting state), where each run is a sequence of mini blocks where subjects perform the VCD-core experiment, + rest periods in between blocks.  Each run approx. 5 minutes
* _Miniblock:_ A series of trials from the same stimulus-task class crossing (one cell in the master table).
* _Trials:_ A series of events where subjects perform the instructed task on the cued stimulus. 
* _Trial events:_ Each trial contains “events”: trial start » spatial cue » ISI » stimulus epoch » (ISI » stimulus epoch) »  response cue » ITI.
* _Trial types:_ We have either single-epoch (one stimulus presentation intervals) or double-epoch (two stimulus presentation intervals, in between is a delay period) trials.
	* _ISI:_ Interstimulus interval (time between two stimulus presentation intervals)
	* _ITI:_ Intertrial interval (time between trials within a miniblock)
	* _IBI:_ Interblock interval (time between mini blocks within a run)

* _Unique images:_ For each stimulus class, we create a unique set of imagees where we manipulate a fix set of features. These manipulations can be at the image level (e.g., gabor contrast level, RDK motion coherence), sometimes they are at the global property level (indoor vs outdoor, manmade vs natural foods). 

* _Semantic Categories:_ For CO and NS we sample different semantic categories:
	* _Superordinate level:_ humans, animals, food, objects, places
	* _Basic level:_ faces, cats, giraffes, tools, houses
	* _Subordinate level:_  *this* cat or *this* dog


## MIT License

Copyright (c) 2024 Eline R. Kupers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
