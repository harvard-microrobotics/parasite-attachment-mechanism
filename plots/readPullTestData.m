function data = readPullTestData(file_path)
    %util function to convert data from csv to array
    data_table = readtable(file_path);
    data = table2array(data_table);
end

