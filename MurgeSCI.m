function [outsci] = MurgeSCI(in_file_name)

    if exist(in_file_name, 'file')

        [path, name, ext] = fileparts(in_file_name);

        post = '_part*.txt';
        final_post = '_Compx.txt';
        outfile_name = strcat(name,final_post);

        findstr = strcat(name, post);

        fileStr= dir(fullfile(path, findstr) );
        files = {fileStr.name}';

        nFile = size(fileStr,1);
        if(nFile > 0)     

            nVertex = size(load(files{1}),1);

            Total_Mat = zeros(nVertex, nFile);

            for i=1:nFile
                Total_Mat(:,i) = load(files{i});
            end

            outsci = sum(Total_Mat,2);

            save(outfile_name,'outsci','-ascii');
        else
            fprintf('No part SCI\n'); 
        end
    else
        fprintf('No OBJ \n');
    end

end