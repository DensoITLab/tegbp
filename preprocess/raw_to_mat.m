% 2022-12-22 jnagata
% convert prophesee raw data to matlab mat data
% require cloning https://github.com/prophesee-ai/prophesee-matlab-scripts
% to outer folder

close all
clearvars

addpath outer/prophesee-matlab-scripts/

%%
rawdata_dir = fullfile('/home/jnagata/petaco/data', 'DNW/221108_EVK2_150mm_f8/500lux');
data_name   = 'recording_2022-11-08_14-02-40';

rawdata_path = fullfile(rawdata_dir, [data_name, '.raw']);
datdata_path = fullfile(rawdata_dir, [data_name, '_cd.dat']);
matdata_path = fullfile(rawdata_dir, [data_name, '_cd.mat']);

if ~exist(datdata_path, 'file')
    command = ['metavision_raw_to_dat -i ', rawdata_path];
    system(command);
    disp('raw to dat : done')
else
    disp('raw to dat : passed')
end

if ~exist(matdata_path, 'file')
    flipX   = 0;
    flipY   = 0;
    event   = load_cd_events(datdata_path, flipX, flipY);
    event.p = logical((event.p+1)/2);
    event.x = uint16(event.x);
    event.y = uint16(event.y);
    save(matdata_path, 'event')
    disp('dat to mat : done')
else
    disp('dat to mat : passed')
end

