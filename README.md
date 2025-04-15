# vcd-stim

MATLAB code repository to create and present the stimuli of the Visual Cognition Dataset (VCD).


## Goal

Create and run the 39 different types of visual experiments of the Visual Cognition Dataset across 30 MRI sessions. Visual experiments are based on crossings of 5 stimulus classes and 9 task classes.

Stimulus classes:
* **Gabors:** Gratings in a gaussian window
* **RDKs:** Random dot motion kinetograms 
* **Dot:** Single simple dot
* **Obj:** Objects
* **NS:** Natural scenes

Task classes:
* **FIX:** Fixation
* **CD:** Contrast-detection
* **SCC:** Superclass categorization
* **PC:** Perceptual categorization
* **WM:** Working memory
* **LTM:** Long-term memory
* **IMG**: Imagery
* **WHAT:** Object categorization
* **WHERE:** Object localization
* **HOW:** Scene/object affordance)

![vcd_master_table](https://github.com/user-attachments/assets/87e1f9ff-ce71-4c62-9548-f9325a09a5c2)


## Dependencies
This stimulus presentation code requires a specific version of [Psychtoolbox-3](https://github.com/Psychtoolbox-3/Psychtoolbox-3) (3.0.14, December 30, 2016), and is written for the iMac that lives permanently at the 7TAS MRI at the Center for Magnetic Resonance Research (CMRR) at the University of Minnesota. 
7TAS iMac specs:
* Retina 5K, 27-inch, 2017
* macOS High Sierra (version 10.13.6), 
* Processor 4.2 GHz Intel Core i7
* Memory 16 GB 2400 MHz DDR4
* Graphics Radeon Pro 575 4096 MB
* MATLAB version 2016b and R2017b
* Psychtoolbox version 3.0.14 (December 30, 2016).

This stimulus presentation code has also been tested on a MacbookPro with macOS 10.14.6 Mojave, and MATLAB version 2018b, PTB version (November 17, 2020), git commit ef093cbf296115badddb995fa06452e34c8c7d02. The stimulus presentation code does not run on MacBook Pro's with Silicon chips (M1 or M2).

## Stimuli and info files

The content of the `workspaces` folder is ignored, as the files are too big. You can view a stored version of the workspace [here](https://drive.google.com/drive/folders/1Boahkioyk5sLrlVFPKiTmLR2RoDlDeaS?usp=sharing).

## Code and folder overview

* `runme_vcdcore.m`	 :	This is the main function to run the core experiment of VCD.
* `vcd_rootPath.m` 	 : 	Function to set the rootpath to relative to the base of this folder.

* `bookkeeping`		 :	Folder with standalone script and vcd related functions to create experimental design and keep track of stimulus conditions during stimulus presentation.
* `external`		 :	Folder with external functions from other toolboxes.
* `params` 			 :	Folder with vcd functions that define display, stimulus, and experimental session parameters.
* `stimpresentation` :	Folder with vcd functions that run stimulus presentation with psychtoolbox.
* `stimcreation`	 : 	Folder with standalone script and functions to create stimuli.
* `tinkers`			 : 	Folder with code tinkering around and such (probably should be removed at some point)
* `utils` 			 : 	Folder with small and simple utility functions.


Folders ignored by git (see .gitignore):
* `workspaces`		:	Folder where the created stimuli, condition tables, and instructions live.
	* `stimuli`		: 	Subfolder where stimuli are stored in matlab (.mat) files
	* `info`		: 	Subfolder where logged stimulus information is stored in csv files, as well as stimulus, display, experiment session params are stored in matlab (.mat) files.	
	* `instructions`: 	Subfolder where instructions are stored in text (.txt) files
* `figs`		    : 	Folder where debug figures are stored
* `data` 		    : 	Folder where subject's button presses and created stimuli are stored (if requested).


## Examples

* Example 1: Create stimuli
  
`s_createStimuli.m`

* Example 2: Create experimental design
  
`s_createDesignMatrix.m`

* Example 3: Present run 01 of session 01 for subject 001:
  
`runme_vcdcore(1,1,1)`

* Example 4: Present run 01 of session 01 for subject 001 in debug mode:
  
`runme_vcdcore(1,1,1, 'debugmode', true)`

* Example 5: Present run 01 of session 01 for subject 001 in debug mode using the psychophysics room monitor:
  
`runme_vcdcore(1,1,1, 'dispName','PPROOM_EIZOFLEXSCAN')`

## Terminology

* `Monitor refresh rate:` The rate by which a monitor runs vertical retrace (VBL). For example, BOLD screen monitor refresh is 120 Hz (for reference: 6 frames = (1000/120) * 6 = 50 ms). Most LCD monitors have a refresh rate of 60 Hz
* `Frames/ framerate / framedur:` Because we present our stimuli at a slower rate than the monitor refresh rate, we “bundle” or “skip” several monitor refreshes into a “frame”. We run our stimuli at 30 Hz, so we bundle 4 120 Hz (BOLDscreen) monitor refreshes or 2 60 Hz (regular LCD) monitor refreshes.
* `Image:`  The retinal image: the image on the BOLD screen that falls onto the retina.
* `Session:` The 2-3 hr scan session at the MRI scanner
* `Run:` Each session consists of approx. 10 runs (+ 2 resting state), where each run is a sequence of mini blocks where subjects perform the VCD-core experiment, + rest periods in between blocks.  Each run approx. 5 minutes
* `Block:` A series of trials from the same stimulus-task class crossing (one cell in the master table).
* `Trials:` A series of events where subjects perform the instructed task on the cued stimulus. 
* `Trial events:` Each trial contains “events”: trial start » spatial cue » ISI » stimulus epoch » (ISI » stimulus epoch) »  response cue » ITI.
* `Trial types:` We have either single-epoch (one stimulus presentation intervals) or double-epoch (two stimulus presentation intervals, in between is a delay period) trials.
	* `ISI:` Interstimulus interval (time between two stimulus presentation intervals)
	* `ITI:` Intertrial interval (time between trials within a miniblock)
	* `IBI:` Interblock interval (time between mini blocks within a run)

* `Unique images:` For each stimulus class, we create a unique set of imagees where we manipulate a fix set of features. These manipulations can be at the image level (e.g., gabor contrast level, RDK motion coherence), sometimes they are at the global property level (indoor vs outdoor, manmade vs natural foods). 

* `Semantic Categories:` For CO and NS we sample different semantic categories:
	* `Superordinate level:` humans, animals, food, objects, places
	* `Basic level:` faces, cats, giraffes, tools, houses
	* `Subordinate level:`  *this* cat or *this* dog

* `Stimulus spatial layout`
  
![vcd_stimulus spatial extend](https://github.com/user-attachments/assets/70146c2c-b355-48f7-afc7-9ec5aefd0851)
Black solid line shows spatial extent of display that is seen by either left or right eye. 
Red dotted line shows spatial extent of display that is seen by both eyes.
Yellow circles are the two peripheral stimulus apertures. Small simple single dot is show at 4.5 deg eccentricity (blue circle)


## MIT License

Copyright (c) 2024 Eline R. Kupers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
