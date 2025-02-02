% https://stackoverflow.com/questions/8981168/running-a-matlab-program-with-arguments
% for now we run exclusively on HK_RR because this is the only dataset supported by the code
% (due to the pre-processed SpeciesNumbering_Standard_HKRR_dataset file loaded below).
% https://stackoverflow.com/questions/3335505/how-can-i-pass-command-line-arguments-to-a-standalone-matlab-executable-running
% https://www.mathworks.com/help/compiler/create-and-install-a-standalone-application-from-matlab-code.html
% https://www.mathworks.com/help/matlab/ref/inputparser.html
% matlab -nodisplay -r "MI_IPA_main(6, 'sub_msa_256.fasta')"
function []=MI_IPA_main(Nincrement, LengthA, replicate, msa_fasta_filename, output_directory)

%This is the main code to run the MI-IPA on the standard HK-RR dataset.
% clear all
close all hidden
set(0,'RecursionLimit',5000)
addpath('Hungarian_algorithm')

%set parameters
% replicate=1;
rng(replicate)
% Nincrement = 6;
% LengthA = 64; %length of first protein (here the HK)

%read data files
% msa_fasta_filename = 'Standard_HKRR_dataset.fasta'; %sequence data file
% msa_fasta_filename = 'sub_msa_256.fasta'
% load SpeciesNumbering_Standard_HKRR_dataset; %read in SpeciesNumbering_extr

%read sequences, adding species number in L+1 and sequence number in L+2
%L is the full length of concatenated sequences, without supplementary indicators such as species and initial index
[encoded_focus_alignment, encoded_focus_alignment_headers, L] = readAlignment_and_NumberSpecies(msa_fasta_filename);
disp(["Concatenated sequence length", L])
%suppress species with one pair
table_count_species =count_species(encoded_focus_alignment);
[encoded_focus_alignment, encoded_focus_alignment_headers] = SuppressSpeciesWithOnePair(encoded_focus_alignment, encoded_focus_alignment_headers, table_count_species);
N = size(encoded_focus_alignment,1); %number of sequences
disp(["Number of sequences with > 1 pair (this may include dummy end and ref start)", N])
%tabulate species and sequences within species
table_count_species =count_species(encoded_focus_alignment);

%number of rounds (last one -> all sequences are in the training set)
Nrounds=ceil(N./Nincrement+1); 
disp(["Number of rounds to perform", Nrounds])

%start from random within-species pairings: scramble the pairings for this.
encoded_training_alignment = ScrambleSeqs(encoded_focus_alignment, LengthA, table_count_species);
%save the species and initial indices of the sequences in the scrambled alignment we start from
filename=strcat(output_directory, '/IniScrambling_Ninc',num2str(Nincrement),'_rep',num2str(replicate),'.txt');
dlmwrite(filename,encoded_training_alignment(:,L+1:L+4),'delimiter','\t')
%in the training set, discard extra indices (species index, initial sequence index)
encoded_training_alignment(:,L+1:L+4)=[];

%initialize
NSeqs_new=0;
Output=zeros(Nrounds,6); 
%Output matrix:
%Each row corresponds to an iteration of the MI-IPA.
%col 1: number of sequences NSeqs in concatenated alignment used as training set
%col 2: effective number of sequences Meff in concatenated alignment used as training set
%col 3: number of TP pairs
%col 4: number of FP pairs
%col 5: number of TP pairs in concatenated alignment used as training set
%col 6: number of FP pairs in concatenated alignment used as training set


%%

% Q: why is there an 'extra' round at the end?
for rounds=1:Nrounds %iterate the process until all sequences are in the training set

    disp(["Round", rounds])
    
    if rounds>1 
        
        %update the training set by adding in the protein pairs with largest gaps predicted at previous round

        %Use the gap to rank pairs
        Results=sortrows(Results,-5);
        
        %number of sequences that will be added to form the training set for this round
        NSeqs_new = NSeqs_new + Nincrement; 
        if NSeqs_new>=size(Results,1)
            NSeqs_new=size(Results,1); %for the last round, all paired sequences will be in the training set
        end
        
        %save to Output the number of TP or FP in the training set
        tps = size(Results(Results(1:NSeqs_new,2)==Results(1:NSeqs_new,3)),1);
        fps = size(Results(Results(1:NSeqs_new,2)~=Results(1:NSeqs_new,3)),1);
        % Question: why is this going wrong...
        disp(["train TPs", tps, "train FPs", fps])
        Output(rounds,5)=tps;
        Output(rounds,6)=fps;

        %construct new training set
        newseqs=zeros(NSeqs_new,L);
        for i=1:NSeqs_new
            newseqs(i,1:LengthA) = encoded_focus_alignment(encoded_focus_alignment(:,L+2)==Results(i,2), 1:LengthA );
            newseqs(i,LengthA+1:L) = encoded_focus_alignment(encoded_focus_alignment(:,L+2)==Results(i,3), LengthA+1:L );  
        end
        encoded_training_alignment = newseqs; 
        
    end

    % things that happen every round
    
    %construct model from training set
    [PMIs, Meff] = Compute_PMIs(encoded_training_alignment,0.15,0.15);
 
    %compute pairings and gap scores for all pairs
    Results =Predict_pairs(encoded_focus_alignment, -PMIs, LengthA, table_count_species);
    % to compute true positives we simply compare the pair index
    tps = size(Results(Results(:,2)==Results(:,3)),1);
    fps = size(Results(Results(:,2)~=Results(:,3)),1);
    disp(["full TPs", tps, "full FPs", fps])
    
    %save the results
    Output(rounds,1)=NSeqs_new;
    Output(rounds,2)=Meff;
    Output(rounds,3)=tps;
    Output(rounds,4)=fps;

end

disp(["TPs", tps, "FPs", fps])


%%
%save Output matrix
filename=strcat(output_directory, '/TP_data_Ninc',num2str(Nincrement),'_rep',num2str(replicate),'.txt');
dlmwrite(filename,Output,'delimiter','\t')

%save the final pairs made and their scores
filename=strcat(output_directory, '/Resf_Ninc',num2str(Nincrement),'_rep',num2str(replicate),'.txt');
dlmwrite(filename,Results,'delimiter','\t')

exit

end