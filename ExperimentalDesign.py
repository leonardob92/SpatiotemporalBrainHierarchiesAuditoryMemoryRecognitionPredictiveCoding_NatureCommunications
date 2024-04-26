
"""

@author: Leonardo Bonetti (with the support of Francesco Carlomagno)

leonardo.bonetti@psych.ox.ac.uk
leonardo.bonetti@clin.au.dk

"""


#################### INFORMATION THAT MUST BE UPDATED BEFORE RUNNING THE SCRIPT ####################

logdir = ('/home/stimuser/Desktop/PatternRecognition_Leonardo/Data_collection_nov_2020') #directory where you have the folder named 'Block_3' which contains two subfolders with the stimuli
frate = 120 #60 #refresh rate of the monitor (i.e.  number of times per second that the screen is redrawn); Typically is 60 Hz in standard monitors and 120 Hz in better monitors. Please check and update according to your setting.
response_buttons = [1, 2] #I assume participants will use buttons 1 and 2 in your response pad. If not, please change them here

#In addition:
  #1) Be sure that the folder 'Block_3' which contains two subfolders with the stimuli is in the 'logdir' that you chose.
  
#################### #################### ####################

#################### Load libraries and set directories ####################

from random import shuffle
from psychopy import prefs
prefs.hardware['audioLib'] = ['PTB']
from psychopy import visual, core, sound, event, gui, logging #, monitors
import os
import numpy as np

### TRIGGER ###
from triggers import setParallelData
setParallelData(0)
### TRIGGER ###

#################### Actual codes ####################

# quit key (to terminate the experiment)
def quitpd():
    win.close()
    logging.flush()
    logfile.close()
    core.quit()

# create folders if not already present
def create_directory(directory_path):
    if not os.path.exists(directory_path):
        try:
            os.makedirs(directory_path)
            print(f"Directory '{directory_path}' created.")
        except OSError as e:
            print(f"Error creating directory '{directory_path}': {e}")
    else:
        print(f"Directory '{directory_path}' already exists.")

event.globalKeys.add(key='escape',func=quitpd)

prd = 1000/frate #conversion screen refresh rate

#GUI for subject ID
ID = gui.Dlg(title = 'subj ID')
ID.addField('ID: ')
subID = ID.show()

#create csv file for log
filename = logdir + '/csv_files/Subj_' + subID[0] + '_Block_3.csv' #name for the csv file
if os.path.exists(filename): #checking that the file does not exist already, to prevent risks of overwriting previous files
    print('THE FILE NAMED ' + 'Subj_' + subID[0] + '_Block_3.csv' + ' ALREADY EXIST!! Check the ID of your current subject')
else:
    #create log_files directory if they do not exist
    directory_path = logdir + '/log_files'
    create_directory(directory_path)
    #create csv_files directory if they do not exist
    directory_path = logdir + '/csv_files'
    create_directory(directory_path)

    logfile = open(filename,'w') #create csv file for the subject
    logfile.write("subj;trial;response;RT;Trig \n") #prepare columns in the csv file
    
    #create log file
    filename = logdir + '/log_files/Subj_' + subID[0] +  '_Block_3.log'
    lastLog = logging.LogFile(filename, level=logging.INFO, filemode='a')
    
    #getting file for learning phase
    os.chdir(logdir + '/Block_3/Learning')
    Bachminor = sound.Sound('Bachminorprelude.wav',volume = 0.5)
    
    #getting files for recognition phase
    os.chdir(logdir + '/Block_3/wav') #changing directory
    pathwave = logdir + '/Block_3/wav' #building path to wave files
    wavefilesd = []
    wavenamesd = []
    for file in sorted(os.listdir(pathwave)): #over files in directory
        if file.endswith(".wav"): #if the file is an audio file
            wavenamesd.append(file) #file name is appended
            wavefilesd.append(sound.Sound(file,volume = 0.5)) #file audio is appended    
    #reproducing several times (here 9 times) the 3 Old excerpts
    wnd = wavenamesd[108:] #extracting Old excerpts (names)
    wnd = wnd * 9 #multiplying them by 9 (names)
    wfd = wavefilesd[108:] #extracting Old excerpts (audio files)
    wfd = wfd * 9 #multiplying them by 9 (audio files)
    wavenames = wavenamesd[0:108] + wnd #getting the final list for the names
    wavefiles = wavefilesd[0:108] + wfd #getting the final list for the files
    wavezip = list(zip(wavenames,wavefiles)) #zipping together files name and files audio
    shuffle(wavezip) #randomising order of trials for the recognition phase
    
    #initializing RT variable to get RTs (not useful in this case, but left it since in variations of this paradigm it would be necessary for quitting purposes in case the def "quitpd" does not work..)
    RT = core.Clock()
    
    #preparing window for the screen, first instruction and fixation cross
    win = visual.Window(fullscr = True, color = 'black') #preparing window
    #learning phase instruction
    instr = visual.TextStim(win,text = 'Learning phase \n\n Now you are going to listen to a complete but short musical piece \n\n Please try to remember it as much as possible \n\n Press 1 to continue',color = 'white')
    fix_c = visual.TextStim(win,text = '+', color = 'white',height = 0.2) #fixation cross

    #actual experiment (learning phase)
    #presenting instruction
    instr.draw()
    win.flip()
    event.waitKeys() #waiting for a key to be pressed
    
    #ACTUAL EXPERIMENT
    #LEARNING PHASE
    for ll in range(2): #over the two repetitions of the musical piece to be encoded
        dumm = 'Repetition ' + str(ll+1) + '\n\n Press 1 to continue' #text for the repetition number 
        playlear = visual.TextStim(win,text = dumm, color = 'white') #preparing information for the message
        playlear.draw() #showing message
        win.flip()
        event.waitKeys() #waiting for a key press whenever the participant is ready
        fix_c.draw() #showing the fixation cross
        win.flip()
        core.wait(0.5) #waiting 0.5 seconds to be sure not to get contamination from the visual presentation of the fixation cross

        ### TRIGGER ###
        nextFlip = win.getFutureFlipTime(clock='ptb') #getting next flip time
        win.callOnFlip(setParallelData, 103) #opening (sending) trigger to the MEG (trigger value = 103)
        Bachminor.play(when = nextFlip) #playing sound, exactly at next flip
        RT.reset()
        resp = None
        for frs in range(int(np.round(50/prd))): #50 ms translated into frames (taking into account the refresh rate of the screen)
            fix_c.draw() #keeping fixation cross
            win.flip()  
        win.callOnFlip(setParallelData, 0) #closing trigger
        for frs in range(int(np.round(50/prd))): #50 ms translated into frames (taking into account the refresh rate of the screen)
            fix_c.draw() #keeping fixation cross
            win.flip()  
        core.wait(24.95) #waiting for the duration of the whole musical piece
    
    #SHORT BREAK
    pausemex = visual.TextStim(win,text = 'Now you have a 20-second break \n\n You can relax :D', color = 'white') #preparing message
    pausemex.draw() #showing message
    win.flip()
    core.wait(20) #waiting 20 seconds
    
    #RECOGNITION PHASE
    #preparing instruction message for recognition phase
    playrec = visual.TextStim(win,text = 'Recognition phase \n\n Now you are going to listen to 135 short musical excerpts \n\n For each of them, please press 1 if the excerpt is "old" (extracted from the musical piece that you have just listened to) or press 2 if it is "new" (not extracted from the musical piece) \n\n Press 1 to continue', color = 'white')
    playrec.draw() #showing message
    win.flip()
    event.waitKeys() #waiting for a key press whenever the participant is ready

    #presentation of stimuli for recognition phase
    for wavve in range(len(wavefiles)): #over melodies for recognition phase
        #displaying the progressive trial number for participants' comfort
        jes = 'trial number ' + str((wavve + 1)) + ' / ' + str(len(wavefiles)) #text for the trial number
        instrrectr = visual.TextStim(win,text = jes,color = 'white') #preparing message
        instrrectr.draw() #showing message
        win.flip()
        core.wait(1) #1 second to allow participants to read the progressive number of the trial
        fix_c.draw() #showing fixation cross
        win.flip()
        core.wait(0.5) #waiting 0.5 seconds to be sure not to get contamination from the visual presentation of the fixation cross

        #preparing trigger value
        if 'old' in wavezip[wavve][0]:  #old melody
            trigval = 10
        else:  #new melody
            trigval = 50
            
        #the following would be a better way to proceed.. here we relied primarily on the .csv file produced by the psychopy code so it was not relevant to further code the melodies.. however, in updated version of this code, we proceeded with the following solution (here commented)
        # #preparing trigger value
        # if 'old' in wavezip[wavve][0]: #old melody
        #     trigval = 100 #trigvalue
        # elif 'new' in wavezip[wavve][0]: #new melody
        #     if 'k1' in wavezip[wavve][0]: #in musical key (weaker variation)
        #         if 't3' in wavezip[wavve][0]: #varied from tone 4
        #             trigval = 110 #trigvalue
        #         elif 't4' in wavezip[wavve][0]: #varied from tone 5
        #             trigval = 120
        #     elif 'k2' in wavezip[wavve][0]: #out of musical key (stronger variation)
        #         if 't3' in wavezip[wavve][0]: #varied from tone 4
        #             trigval = 130
        #         elif 't4' in wavezip[wavve][0]: #varied from tone 5
        #             trigval = 140    
            
        ### TRIGGER ###    
        nextFlip = win.getFutureFlipTime(clock='ptb') #getting next flip time
        win.callOnFlip(setParallelData, trigval) #opening (sending) trigger to MEG
        ### TRIGGER ###
        
        wavezip[wavve][1].play(when = nextFlip) #playing sound, exactly at next flip
        event.clearEvents(eventType='keyboard') #making sure that no keys are recorded before the beginning of the musical sequence (in case, e.g., the participant presses a key randomly between different trials)
        RT.reset()
        resp = None
        
        ### TRIGGER ###
        for frs in range(int(np.round(50/prd))): #time frames corresponding to 50 ms
            fix_c.draw()
            win.flip() 
        win.callOnFlip(setParallelData, 0) #closing trigger
        for frs in range(int(np.round(50/prd))): #time frames corresponding to 50 ms
            fix_c.draw()
            win.flip()
        ### TRIGGER ###

        while resp == None: #while there is no response
            key = event.getKeys(keyList = [str(response_buttons[0]),str(response_buttons[1]),'escape']) #looking for response (either one of the two keys used for recording participant's response or the escape key in case you need to quit the experiment in case the def "quitpd" does not work)
            if len(key) > 0: #if response is given
                rt = RT.getTime() #getting RT
                resp = key[0][0] #getting the actual response (e.g. '1' or '2')
            elif RT.getTime() > 3.7: #otherwise waiting 3.7 seconds - maximum waiting time if participant does not reply
                resp = 0 #in this case response is coded as '0'
                rt = RT.getTime() #getting fixed (maximum) RT
        for frs in range(int(np.round(1950/prd))): #time frames corresponding to 1950 ms
            fix_c.draw()
            win.flip()
        logging.flush()           
        #writing SUBJ ID, trial ID, subject's response, RT, trigger value
        lrow = '{};{};{};{};{} \n' #new row in the csv file
        lrow = lrow.format(subID[0],wavezip[wavve][0], resp, round(rt*1000),str(trigval)) #storing behavioural information
        logfile.write(lrow) #writing on csv
        if 'escape' in key: #just in case experimenter presses 'escape' on keyboard because they want to quit the experiment (and the def "quitpd" does not work)
            logging.flush()
            core.quit()
    
    #closing csv file
    logfile.close()   
    
    #final message
    playlear = visual.TextStim(win,text = 'Thank you very much! \n\n The next experimental block will start soon', color = 'white')
    playlear.draw()
    win.flip()
    event.waitKeys(keyList = 'space')

####################
