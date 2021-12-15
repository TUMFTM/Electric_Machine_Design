% -------------------------------------------------------------------------
% TU Munich - Institute of Automotive Technology
% -------------------------------------------------------------------------
% Modell for the design and analysis of PMSM or ASM (MEAPA)
% -------------------------------------------------------------------------
% Autor:    Svenja Kalt (svenja.kalt@tum.de), 
%           Jonathan Erhard 
% -------------------------------------------------------------------------

%% NOTE:
% For the update of the headers in all m-files all m-files of the
% MEAPA tool must be closed. Then in this script under
% newString the new version with editor code is entered.
% Afterwards the m-file 'UpdateVersion' is saved and also % closed.
% closed. Finally the command 'UpdateVersion' must be % entered in the command window.
% (without quotation marks) must be executed.

%% Parameters
% The directory in which to replace files. Currently this code does not modify files in
% sub-directories
directory           = pwd;
% The string that will be replaced
oldString           = sprintf('v0.8 - 13.06.2019 - je');
% The replacement string
newString           = sprintf('v0.8 - 13.06.2019 - je');
% The file name condition - what type of files will be examined
% It must contain any of the English character set (letters, numbers or underscore
% character i.e. a-zA-Z_0-9) and ends with a ".m" MATLAB extension (use \.txt for text files)
regularExpression   = '[\w]+\.m';

%% Determine files to update, and update them as necessary
% Change the current directory to the user-specified one
cd(directory);
% Put the details of all files and folders in that current directory into a structure
allFilesInDirectory = dir;
allFilesInDirectory = allFilesInDirectory(~ismember({allFilesInDirectory.name},{'.','..'}));
% Initialise indexes for files that do and do not contain oldString
filesWithStringIndex = 1;
filesWithoutStringIndex = 1;
filesWithString{filesWithStringIndex} = '';
filesWithoutString{filesWithoutStringIndex} = '';
% For the number of files and folders in the directory
for idx = 1 : length(allFilesInDirectory)

    % If the file name contains any of the English character set (letters, numbers or
    % underscore character i.e. a-zA-Z_0-9) and ends with a ".m" filetype...
    if (~isempty ( regexp(allFilesInDirectory(idx).name, '[\w]+\.m','match') ))

        % Open the file for reading
        fileIdRead  = fopen([allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name], 'r');

        % Extract the text
        fileText = fscanf(fileIdRead,'%c');

        % Close the file
        fclose(fileIdRead);

        % Search for occurrences of oldString
        occurrences = strfind(fileText,oldString);

        % If an occurrence is found...
        if ~isempty(occurrences)

            % Replace any occurrences of oldString with newString
            fileTextNew = strrep(fileText, oldString, newString);

            % Open the file for writing
            fileIdWrite = fopen([allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name], 'w');

            % Write the modified text
            fprintf(fileIdWrite, '%c', fileTextNew);

            % Close the file
            fclose(fileIdWrite);

            % Update the list of files that contained oldString
            filesWithString{filesWithStringIndex} = allFilesInDirectory(idx).name;

            % Update the index for files that contained oldString
            filesWithStringIndex = filesWithStringIndex + 1;

        else
            % Update the list of files that did not contain oldString
            filesWithoutString{filesWithoutStringIndex} = allFilesInDirectory(idx).name;

            % Update the index for files that did not contain oldString
            filesWithoutStringIndex = filesWithoutStringIndex + 1;

        end
    elseif(allFilesInDirectory(idx).isdir)
        Folder = [allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name];
        [filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex] = Recursive(Folder,filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex,newString,oldString);
    end
end
%% Display what files were changed, and what were not
% If the variable filesWithString exists in the workspace
if exist('filesWithString','var')
    disp('Aktualisierte Dateien:');
    % Display their names
    for i = 1:filesWithStringIndex-1, disp(filesWithString{i}); end
% else
%     disp('No files contained the target string');
end
% Insert a clear line between lists
% disp(' ');
% If the variable fileWithoutString exists in the workspace
% if exist('filesWithoutString','var')
%     disp('Files that were not updated:');
%     % Display their names
%     for j = 1:filesWithoutStringIndex-1, disp(filesWithoutString{j}); end
% else
%     disp('All files contained the target string.');
% end

clearvars ans filesWithString filesWithStringIndex filesWithoutString filesWithoutStringIndex Folder allFilesInDirectory idx fileIdWrite fileIdRead fileTextNew newString oldString directory regularExpression occurrences fileText i j

%% Recursive function for subfolders
function [filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex] = Recursive(Folder,filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex,newString,oldString)
    allFilesInDirectory = dir(Folder);
    allFilesInDirectory = allFilesInDirectory(~ismember({allFilesInDirectory.name},{'.','..'}));
    
    for idx = 1 : length(allFilesInDirectory)

        % If the file name contains any of the English character set (letters, numbers or
        % underscore character i.e. a-zA-Z_0-9) and ends with a ".m" filetype...
        if (~isempty ( regexp(allFilesInDirectory(idx).name, '[\w]+\.m','match') ))

            % Open the file for reading
            fileIdRead  = fopen([allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name], 'r');

            % Extract the text
            fileText = fscanf(fileIdRead,'%c');

            % Close the file
            fclose(fileIdRead);

            % Search for occurrences of oldString
            occurrences = strfind(fileText,oldString);

            % If an occurrence is found...
            if ~isempty(occurrences)

                % Replace any occurrences of oldString with newString
                fileTextNew = strrep(fileText, oldString, newString);

                % Open the file for writing
                fileIdWrite = fopen([allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name], 'w');

                % Write the modified text
                fprintf(fileIdWrite, '%c', fileTextNew);

                % Close the file
                fclose(fileIdWrite);

                % Update the list of files that contained oldString
                filesWithString{filesWithStringIndex} = allFilesInDirectory(idx).name;

                % Update the index for files that contained oldString
                filesWithStringIndex = filesWithStringIndex + 1;

            else
                % Update the list of files that did not contain oldString
                filesWithoutString{filesWithoutStringIndex} = allFilesInDirectory(idx).name;

                % Update the index for files that did not contain oldString
                filesWithoutStringIndex = filesWithoutStringIndex + 1;

            end
        elseif(allFilesInDirectory(idx).isdir)
            Folder = [allFilesInDirectory(idx).folder '/' allFilesInDirectory(idx).name];
            [filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex] = Recursive(Folder,filesWithString,filesWithStringIndex,filesWithoutString,filesWithoutStringIndex,newString,oldString);
        end
    end
end