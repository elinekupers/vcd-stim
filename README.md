# vcd-stim

MATLAB code repository to create and present the stimuli of the Visual Cognition Dataset (VCD).


## Goal

Create and run the 32 different types of visual experiments of the Visual Cognition Dataset across 27 MRI sessions. Visual experiments are based on crossings of 5 stimulus classes and 9 task classes.

Stimulus classes:
* **Gabors:** Gratings in a gaussian window
* **RDKs:** Random dot motion kinetograms 
* **Dot:** Single simple dot
* **Obj:** Objects
* **NS:** Natural scenes

Task classes:
* **FIX:** Fixation brightness task 
* **CD:** Contrast change task
* **SCC:** Superclass Categorization task
* **PC:** Perceptual categorization:
	* Gabor tilt task
	* RDK motion direction task
	* Dot position task
	* Object rotation task
	* Scene Indoor/outdoor task
* **WM:** Working memory
	* Gabor tilt memory task
 	* RDK motion direction memory task
  	* Dot position memory task
  	* Object rotation memory task
  	* Scene memory task 
* **LTM:** Long-term memory matching task
* **IMG**: Imagery task
* **WHAT:** Object categorization task
* **WHERE:** Object localization task
* **HOW:** Scene/object affordance task

![vcd_master_table](https://github.com/user-attachments/assets/c50f7032-c08f-4271-ad5b-e5fa9a11a89d)



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
Base functions:
* `vcd_rootPath.m` : Function to set the rootpath to relative to the base of this folder.
* `vcd.m` : Function to get basic info about the stimulus classes, task classes, and stimuli used in VCD.
* `runme_vcdcore.m` : This is the main function to run the core experiment of VCD.

	The `runme_vcdcore.m` mandatory inputs are (in this order): 
		* subject number : 1-999 (nr will be zero padded)
		* session number : 1-27 for MRI; 1 for behavior.
		* session type   : 1 for version A; 2 for version B. Version B only exists for MRI session 1 and 27).
		* run_nr 		 : 1-10 for MRI; 1-13 for behavior.
		* display name   :'7TAS_BOLDSCREEN32' for BOLDscreen monitor in the 7TAS; 'PPROOM_EIZOFLEXSCAN' for Eizo Flexscan monitor in the psychophysics room).
 
	The `runme_vcdcore.m` can take several optional inputs, here are the most important ones:
		* randomization_file : path to file where subject's time_table_master lives. If randomization_file is left empty, the stimulus presentation code will generate a new file from the corresponding condition_master mat file for that experiment/display. Example path: `fullfile(vcd_rootPath', 'data', 'BEHAVIOR', 'vcd_subj000', 'vcd_subj000_time_table_master_PPROOM_EIZOFLEXSCAN_20250609T200011.mat)` 
		* `wantsynctest`     : Logical flag to run Psychtoolbox sync test or not. Default: true.
		* `wanteyetracking`  : Logical flag to initialize Eyelink and record eyetracking data. Default: false.
		* `ptbMaxVBLstd`     : (double) time in seconds referring to the max standard devation across VBL time stamps collected during the PTB sync test. If standard deviation happens to be larger than this value, the sync test fails. Default: 0.0006 s (or 0.6 ms).

 	The `runme_vcdcore.m` outputs are (in this order): 
 		* data  	 : (struct) behavioral data including button presses and their time stamps, VBL time stamps, subject's run table and frames describing every event during the experimental run. We also store other variables present in the workspace (except for uint8 images).
 		* all_images : (struct) all uint8 stimuli and alpha transparency masks used in the vcd-core experiment. One may want to have all stimuli ready in the workspace to avoid reloading the stimulus mat-files every single run.


Folders:
* `analysis` : Folder with vcd related functions to analyze behavioral and eye gaze data.
* `bookkeeping` : Folder with standalone script and vcd related functions to create experimental design and keep track of stimulus conditions during stimulus presentation.
* `external` : Folder with external functions from other toolboxes.
* `params` : Folder with vcd functions that define display, stimulus, and experimental session parameters.
* `stimpresentation` : Folder with vcd functions that run stimulus presentation with psychtoolbox.
* `stimcreation` : Folder with standalone script and functions to create stimuli.
* `tinkers` : Folder with code tinkering around and such (probably should be removed at some point)
* `utils` : Folder with small and simple utility functions.


Folders with contents that are mostly ignored by git (see .gitignore):
* `workspaces` : Folder where the created stimuli, condition tables, and instructions live.
	* `stimuli` : Subfolder where stimuli are stored in matlab (.mat) files
	* `info` : Subfolder where logged stimulus information is stored in csv files, as well as stimulus, display, experiment session params are stored in matlab (.mat) files. (TWO EXCEPTIONS: we do not ignore the condition_master_*.mat files as this file is used by vcd.m, and we do not ignore the .csv info files as they are small and contain important information about the stimuli.)
	* `instructions` : Subfolder where instructions are stored in text (.txt) files
* `figs` : Folder where debug figures are stored
* `data` : Folder where subject's button presses, created run_table and run_frames are stored.


## Examples

Example 0: Get all stimulus class names

`vcd('stimulusclassnames',[])`


Example 1: Create stimuli

`s_createStimuli.m`



Example 2: Create experimental design

`s_createDesignMatrix.m`


Example 3: For subject 001, present MRI session 01 (wide), session type A, run 01, with eye tracking:
  
`runme_vcdcore(1,1,1,1,'7TAS_BOLDSCREEN32','wanteyetracking', true)`


Example 4: For subject 001, present MRI session 01 (wide), session type B, run 01, with eye tracking:
  
`runme_vcdcore(1,1,2,1,'7TAS_BOLDSCREEN32','wanteyetracking', true)`


Example 5: For subject 001, present behavioral session 01, session type A, run 02, without eye tracking (note: there is only type A for behavior):
  
`runme_vcdcore(1,1,1,2,'PPROOM_EIZOFLEXSCAN')`


Example 6: For subject 002, present MRI session 07 (deep), session type A, run 01, without eye tracking (note: there is only type A for MRI session 2-26):
  
`runme_vcdcore(2,7,1,5,'7TAS_BOLDSCREEN32')`



## Terminology

* `Monitor refresh rate:` The rate by which a monitor runs vertical retrace (VBL). For example, BOLD screen monitor refresh is (approximately) 120 Hz (for reference: 6 frames = (1000/120) x 6 = 50 ms). Most standard monitors have a refresh rate of about 60 Hz. But other rates are possible, like 75 Hz, 85 Hz.
* `Frames / frame rate / frame duration:` Because we present our stimuli at a slower rate than the monitor refresh rate, we “bundle” or “skip” several monitor refreshes into a “frame”. We run our stimuli at 30 Hz, so we bundle 4 120 Hz (BOLDscreen) monitor refreshes or 2 60 Hz (regular LCD) monitor refreshes.
* `Stimulus:` Either the thing on the left, right, or centrally positioned. There can be two stimuli on the screen simultaneously.
* `Session:` The 2-3 hr scan session at the MRI scanner
* `Run:` Each session consists of approx. 10 runs (+ 2 resting state), where each run is a sequence of blocks where subjects perform the VCD-core experiment, + rest periods in between blocks.  Each run approx. 5 minutes. Each run also includes an eyetracking block at the beginning.
* `Block:` A series of trials from the same stimulus-task class crossing (one cell in the master table).
* `Trials:` A series of events where subjects perform the instructed task on the cued stimulus. 
* `Trial events:` Each trial contains “events”: trial start » spatial cue » ISI » stimulus epoch » (ISI » stimulus epoch) » response cue » ITI.
* `Fixation circle:` it gets brighter or dimmer.
* `Fixation circle rim:` Part of the fixation circle. The outer rim is what indicates trial starting and ending. It can get thick or thin.
* `Spatial cue:` Either left, right, or central (using red to indicate).
	* `Cued location:` either left, right, or central.
	* `Uncued location:` typically refers to right/left if the cued location is left/right.
* `Task cue:` A short screen of text at the beginning of every block.
* `Trial types:` We have either one-image trials (one stimulus presentation, the stimulus) or two-image trials (two stimulus presentation (WM, LTM: reference image -> test image; IMG: text prompt -> test image), in between is a delay period).
* `Delay period:` For two-image trials, the time between the offset of the first stimulus presentation and the onset of the second stimulus presentation
* `ITI:` Intertrial interval (time between the offset of the response window (rim becomes thin) and the onset of the rim thickening of the next trial within the miniblock)
* `IBI:` Interblock interval (time between the offset of the response window of the last trial in a block, and the onset of the task instruction window of the following miniblock). IBIs involve just the fixation circle.
* `Correct image:` The image that you should have learned to associate with the reference image) vs. incorrect image (some incorrect images are lures (applies only to LTM; defined as incorrect associated images that are within stimulus class (in the classic group) or incorrect similarly-looking images (in the natural group)), modified stimulus (applies only NS-WM), stimulus pair (AB) where B is associated with A (LTM), test dots (IMG).
* `Eyetracking terms:` targets, eyetracking block, Eyelink calibration
* `Rest time:` Whenever there is “nothing” for the participant to do other than fixate.
* `Stimulus-task crossing:` the crossing between a stimulus class and task class (e.g., NS-WM).
* `Unique stimuli:` All the images/movies within a single stimulus class. The number of unique stimuli is determined by the number of fully-crossed manipulated stimulus features (up to 3 stimulus features, see below). The manipulated stimulus features of interest can be at the image level (e.g., gabor contrast level, RDK motion coherence), or at the global property level (indoor vs outdoor, manmade vs natural foods). 
* `Semantic Categories:` For CO and NS we sample different semantic categories:
	* `Superordinate level:` humans, animals, food, objects, places
	* `Basic level:` faces, cats, giraffes, tools, houses
	* `Subordinate level:`  *this* cat or *this* dog

## Stimulus spatial layout
* `Background` exists of a full-field pink noise image and a mean-luminance gray "puzzle piece" cut out. The cutout is the union of all types of stimuli plus a 1-degree buffer zone on each side.
* `Gabors, RDKs, and Objects` are presented within circular apertures (4 degree diameter), centered at 4 degrees eccentricity on the horizontal meridian. We choose these parafoveal locations to ensure subjects can see both left and right stimuli with both eyes.
* `Single dots` are each 1 degree diameter and live on an iso-eccentric ring at 4 degrees eccentricity (to match the Gabors, RDKs, Objects).
* `Scenes` are 8.4 degree width and height, and are presented centrally.
* `Fixation circle` consists of two parts: an inner part (the part that changes luminance, 0.14 degrees in diameter) and a rim. During ITIs or IBIs, the rim is thin (~0.20 degrees diameter). During a trial, the rim is thick (~0.25 degrees diameter). The fixation circle is presented on top the background (or natural scene) and is 50% transparent.
* `Alpha transparency masks` are applied to each element of the display to obscure the edges of each stimulus' support (e.g., a circular gabor lives on a rectangular support). One exception is scenes, which do not require an alpha mask.
  	  
![Gabor stimulus display for 7TAS BOLDScreen (same for RDK)](https://github.com/user-attachments/assets/dd822b0b-ec01-49e5-b27e-a745d9b36ef7)

![Object stimulus display for 7TAS BOLDScreen](https://github.com/user-attachments/assets/5b5dcec0-75da-44fc-b3f4-8fca9e55435a)

![Single dot stimulus display for 7TAS BOLDScreen](https://github.com/user-attachments/assets/7d5c8672-d7d1-4f4c-a630-68f52c8ced24)

![Natural scenes stimulus display for 7TAS BOLDScreen](https://github.com/user-attachments/assets/bcdbf59c-0b6c-4e26-b739-d765a6874aca)




## MIT License

Copyright (c) 2025 Eline R. Kupers

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
