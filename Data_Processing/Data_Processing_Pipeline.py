
from Data_Processing.EAG_DataProcessing_Library import *

#designate which folders need to be analyzed. this should be the experiment level directory
#Dir=['Manduca2Reciever_Data/Raw_Data/WT/062723/062723M1fA1_WT', 'Manduca2Reciever_Data/Raw_Data/WT/070623/070623M1fA2_WT'
   #,'Manduca2Reciever_Data/Raw_Data/WT/070623/070623M2fA2_WT']

Dir=['Manduca2Reciever_Data/Raw_Data/OR6KO/062723/062723M1fA1_OR6KO',
     'Manduca2Reciever_Data/Raw_Data/OR6KO/062823/OR6KO/062823M1fA2_OR6KO','Manduca2Reciever_Data/Raw_Data/OR6KO/062823/OR6KO/062823M2fA2_OR6KO',
    'Manduca2Reciever_Data/Raw_Data/WT/062723/062723M1fA1_WT', 'Manduca2Reciever_Data/Raw_Data/WT/070623/070623M1fA2_WT'
    , 'Manduca2Reciever_Data/Raw_Data/WT/070623/070623M2fA2_WT']

process_data(DL=Dir, record_channels=['EAG1','EAG2'], savedir='Manduca2Reciever_Data/Buttered/Norm_Independent/',
             SUM=False, Norm='YY', Smoothen=False, LOG=False, Butter=[1,5,2], B_filt=True, RETURN='SAVE')

process_data(DL=Dir, record_channels=['EAG1','EAG2'], savedir='Manduca2Reciever_Data/Buttered/Norm_Sum/',
             SUM=True, Norm='YY', Smoothen=False, LOG=False, Butter=[1,5,2], B_filt=True, RETURN='SAVE')

process_data(DL=Dir, record_channels=['EAG1','EAG2'], savedir='Manduca2Reciever_Data/Buttered/Independent/',
             SUM=False, Norm=False, Smoothen=False, LOG=False, Butter=[1,5,2], B_filt=True, RETURN='SAVE')

process_data(DL=Dir, record_channels=['EAG1','EAG2'], savedir='Manduca2Reciever_Data/Buttered/Sum/',
             SUM=True, Norm=False, Smoothen=False, LOG=False, Butter=[1,5,2], B_filt=True, RETURN='SAVE')


