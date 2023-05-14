function myaudiowrite(filename,y,Fs,varargin)
%AUDIOWRITE write audio files
%   AUDIOWRITE(FILENAME,Y,FS)  writes data Y to an audio
%   file specified by the file name FILENAME, with a sample rate
%   of FS Hz.
%   Stereo data should be specified as a matrix with two columns.
%   Multi-channel data should be specified as a matrix of N columns.
%
%   AUDIOWRITE(URL, Y, FS) writes the audio file at a remote location
%   (URL). When writing data to remote locations, you must specify the full
%   path using a uniform resource locator (URL). For example,
%   to write an audio file to Amazon S3 cloud specify the full URL for the
%   file:
%       s3://bucketname/path_to_file/my_audio.wav
%
%   For more information on accessing remote data, see "Work with Remote Data"
%   in the documentation.
%   The file format to use when writing the file is inferred from
%   FILENAME's extension.  Supported formats are as follows:
%
%   Format         File Extension(s)   Compression Method 
%   ---------------------------------------------------------
%   Wave           .wav                None                
%   MPEG-4 Audio   .m4a, .mp4          AAC                 
%   FLAC           .flac               FLAC (Lossless)     
%   Ogg/Vorbis     .ogg, .oga          Vorbis          
%   Ogg/Opus       .opus               Opus
%   MP3            .mp3                MPEG-1/2/2.5 LAYER-III
%
%   AUDIOWRITE(FILENAME,Y,FS,Name1,Value1,Name2,Value2, ...) 
%   Specifies optional comma-separated pairs of Name,Value arguments, 
%   where Name is the argument name and Value is the corresponding value. 
%   Name must appear inside single quotes (' '). You can specify several 
%   name and value pair arguments in any order as Name1,Value1,...,NameN,
%   ValueN.  Valid Name,Value arguments are as follows:
%
%   'BitsPerSample'  Number of bits per sample to write out to the audio 
%                    file.  
%                    Only supported for WAVE (.wav) and FLAC(.flac) files.
%                    Valid values are 8,16,24,32, or 64 and vary depending 
%                    upon the supported format. 
%                    See "Output Data Type" section below for more details.
%
%   'BitRate'        Number of kilobits per second (kbps) used for compressed 
%                    audio files.  In general, the larger the BitRate The 
%                    higher the compressed audio quality.
%                    Only Supported for MPEG-4 Audio (.m4a, .mp4) files.
%
%   'Quality'        A scalar integer in the range [0,100], where 0 is lower quality
%                    and higher compression, and 100 is higher quality and 
%                    lower compression. Only available for OGG (.ogg),
%                    OPUS (.opus), and MP3 (.mp3).
%
%   'Title'          Character vector or string scalar representing a title
%                    to be written to the file.
%
%   'Artist'         Character vector or string scalar representing the
%                    artist or author to be written to the file.
%
%   'Comment'        Character vector or string scalar representing a
%                    comment to be written to the file.
%
%   Input Data Ranges  
%   Y is specified as an m-by-n matrix, where m is the
%   number of audio samples to write and n is the number of audio channels 
%   to write.  If either m or n is 1, then audiowrite assumes that 
%   dimension specifies the number of audio channels, and the other 
%   dimension specifies the number of audio samples.
%   
%   The valid range for the data in Y depends on the data type of Y.
%   Supported data types are as follows:
% 
%   Data Type of Y    Valid Range for Y
%   -----------------------------------
%       uint8            0 <= Y <= 255
%       int16       -32768 <= Y <= +32767
%       int32        -2^31 <= Y <= 2^31-1
%       single        -1.0 <= Y <= +1.0
%       double        -1.0 <= Y <= +1.0
%
%   Data beyond the valid range is clipped.
% 
%   If Y is single or double, then audio data in Y should be normalized 
%   to values in the range -1.0 and 1.0, inclusive.
%
%   Output Data Type
%   The native data type written to the audio file is determined by 
%   the file format, the data type of Y, and the specified output 
%   BitsPerSample.
% 
%   File Formats  Data Type of Y     Output BitsPerSample  Output Data Type
%   -----------------------------------------------------------------------
%   WAVE (.wav)	  uint8,int16,int32         8                  uint8
%                 single,double             16                 int16
%                                           24                 int32
%                 ---------------------------------------------------------
%                 uint8,int16,int32         32                 int32
%                 ---------------------------------------------------------
%                 single,double             32                 single
%                 ---------------------------------------------------------
%                 single,double             64                 double
%   -----------------------------------------------------------------------
%   FLAC (.flac)  uint8,int16,int32,        8                  int8
%                 single,double             16                 int16
%                                           24                 int32
%   -----------------------------------------------------------------------
%   MPEG-4        uint8,int16,int32,        N/A                single
%   (.m4a,.mp4),  single,double
%   OGG (.ogg,
%   .oga,.opus),
%   MP3 (.mp3)	
%   -----------------------------------------------------------------------
%
% Example
%     % Write audio to a local file
%     audiowrite('Local_folder/sample_audio.wav', y, Fs);
%     
%     % Write audio to an Amazon S3 location
%     audiowrite('s3://bucketname/path_to_file/sample_audio.wav', y, Fs);
%     
%   See also AUDIOINFO, AUDIOREAD

%   Copyright 2012-2022 The MathWorks, Inc.


if nargin > 0
    filename = convertStringsToChars(filename);
end

if nargin > 3
    [varargin{:}] = convertStringsToChars(varargin{:});
end

narginchk(3,13);

% Parse inputs:
props = parseInputs(filename, y, Fs, varargin);

% Check if filename is a URL or not
validURL = matlab.io.internal.vfs.validators.isIRI(filename);
validPrefix = matlab.io.internal.vfs.validators.hasIriPrefix(filename);

if validURL && validPrefix
   % check whether credentials are set
    import matlab.io.internal.vfs.validators.validateCloudEnvVariables;
    validateCloudEnvVariables(filename);
    
    % Get the correct filename
    filenameWithoutPath = matlab.io.internal.vfs.validators.IRIFilename(filename);
    remoteFolder = extractBefore(filename, filenameWithoutPath);
    [~,~, ext] = fileparts(filenameWithoutPath);
    if strncmp(ext,'.',1)
        ext = ext(2:end);
    end
    
    validateExtension(ext, filename);
    
    localToRemote = matlab.io.internal.vfs.stream.LocalToRemote(remoteFolder);
    % Set the file name that needs to go into the remote folder.
    setfilename(localToRemote,filenameWithoutPath);
    props.filename = localToRemote.CurrentLocalFilePath;
    
    
    validateRemotePath(localToRemote);
    % Setting the bool variable 'fileExisted' as 1 so that the temp file is
    % not deleted in 'handleException', and the temp file is always cleaned
    % by the vfs cleanup function
    fileExisted = 1;
     
elseif validPrefix && ~validURL
    % possibly something wrong with the URL such as a malformed Azure
	% URL missing the container@account part
    error(message("MATLAB:virtualfileio:stream:invalidFilename", filename));
else
    [props.filename, fileExisted] = validateFilename( props.filename );
end


% audiowrite supports real values only
 validateattributes(y,{'numeric'},{'real'});

% y can only be 1D or 2D
if ~ismatrix(y) 
  error(message('MATLAB:audiovideo:audiowrite:invalidDimensions'));
end

% If input is a vector, force it to be a column:
if size(y,1)==1
   y = y(:);
end

[~, props.Channels] = size(y);

import multimedia.internal.audio.file.PluginManager;

try
    writePlugin = PluginManager.getInstance.getPluginForWrite(props.filename);
catch exception
    % The exception has been fully formed. Only the prefix has to be
    % replaced.
    exception = PluginManager.replacePluginExceptionPrefix(exception, 'MATLAB:audiovideo:audiowrite');
    
    throw(exception);
end
    
import multimedia.internal.audio.file.FrameSize;
try
    options.Filename = props.filename;
    
    % Create Channel object
    channel = matlabshared.asyncio.internal.Channel( ...
        writePlugin,...
        PluginManager.getInstance.MLConverter, ...
        Options = options, ...
        StreamLimits = [0, 0]);
    
    % The channel exposes settable custom properties depending upon the 
    % file type.  Validate that the input properties are valid for this
    % type of file.
    props = validateProperty('BitRate',props,'BitRate',channel,128);
    props = validateProperty('BitsPerSample',props,'FileDataType',channel,16);    
    props = validateProperty('Quality',props,'Quality',channel,75);  
    
    
    % create and configure a transform filter to convert 
    % the incoming data into a form that the channel's device plugin
    % will expect.
    channel.OutputStream.addFilter( ...
        PluginManager.getInstance.TransformFilter, ...
        []);
    
    options.FilterOutputDataType = filterOutputDataTypeFromFile(...
        props.filename, ...
        class(y), ...
        channel.DataType);
    
    options.FilterTransformType = 'InterleaveTranspose';
    options.SampleRate = double(props.Fs);
    options.NumberOfChannels = double(props.Channels);
    options.BitRate = double(props.BitRate) * 1000; % convert to Bits-Per-Second
    options.Quality = double(props.Quality) / 100; % value must be between 0 and 1
    options.FileDataType = fileDataTypeFromBitsPerSample(...
        props.BitsPerSample,...
        class(y),...
        props.filename);
    
    options.Comment = props.Comment;
    options.Artist = props.Artist;
    options.Title = props.Title;
    
    % Just prior to writing y to the input range if needed
    %y = clipInputData(y, options.FileDataType);
    
    channel.open(options);
 
    channel.OutputStream.write(y, double(FrameSize.Optimal));
    
    % Warn about writing Mp4 metadata on Mac
    warnIfMp4Metadata(props);
    
    channel.close();
catch exception
    try
        channel.close(); % Close the channel
    catch ME
        handleException(ME, props, fileExisted);
    end
    handleException(exception, props, fileExisted);
end

if validURL
    upload(localToRemote);
end
end

function props = parseInputs(filename, y, Fs, pvpairs)

p = inputParser;
p.addRequired('filename',@(x)validateattributes(x,{'char','string'},{'nonempty'}));

p.addRequired('y',...
    @(x)validateattributes( x, ...
    {'uint8','int16','int32','single','double'},{'nonempty'}));

p.addRequired('Fs',@(x)validateattributes(x,{'numeric'},...
    {'nonempty','positive','integer','nonnan'}));

p.addParameter('BitsPerSample',[],@(x)validateattributes(x,{'numeric'},...
    {'nonempty','positive','integer'}));

p.addParameter('BitRate',[],@(x)validateattributes(x,{'numeric'},...
    {'nonempty','positive','integer'}));

p.addParameter('Quality',[],@(x)validateattributes(x,{'numeric'},...
    {'nonempty','integer','>=',0,'<=',100}));

p.addParameter('Title',[],@(x)validateattributes(x,{'char','string'},{}));
p.addParameter('Artist',[],@(x)validateattributes(x,{'char','string'},{}));
p.addParameter('Comment',[],@(x)validateattributes(x,{'char','string'},{}));

p.CaseSensitive = false;
p.KeepUnmatched = true;
p.FunctionName='audiowrite';

parse(p,filename, y, Fs, pvpairs{:});

% Partially Match 'unmatched' properties
if ~isempty(p.Unmatched)
    validParams = { ...
        'BitsPerSample',...
        'BitRate',...
        'Title',...
        'Artist',...
        'Comment'};
    
    props = p.Results;
    
    unmatchedParams = fieldnames(p.Unmatched);
    for ii=1:length(unmatchedParams)
        matchedParam = validatestring(unmatchedParams{ii}, validParams,'audiowrite');
        
        props.(matchedParam) = p.Unmatched.(unmatchedParams{ii});
    end
end


end

function dataType = fileDataTypeFromBitsPerSample(bitsPerSample, inputDataType, filename)

validBitsPerSample = [8 16 24 32 64];
validDataTypes = {'uint8' 'int16' 'int24' 'single' 'double'};

% Flac files only support 8, 16, and 24 bit
filepath = multimedia.internal.io.FilePath(filename);
if strcmpi(filepath.Extension, 'flac')
    validBitsPerSample = validBitsPerSample(1:3);
    validDataTypes = validDataTypes(1:3);
end

if ~ismember(bitsPerSample,validBitsPerSample)
    values = sprintf('%d,',validBitsPerSample);
    values = values(1:end-1); % remove trailing ','
    error(message(...
        'MATLAB:audiovideo:audiowrite:invalidBitsPerSample',...
        sprintf('%d',bitsPerSample), values));
end

dataType = validDataTypes{bitsPerSample == validBitsPerSample};

if bitsPerSample == 32 && ismember(inputDataType,{'uint8','int16','int32'})
    % If bitsPerSample is 32, Integer input data should be 
    % int32 instead of single.
    dataType = 'int32';
end
end

function dataType = filterOutputDataTypeFromFile(filename, inputDataType, channelDataType)

% For most files, just use the underlying Channel's device plugin DataType
dataType = channelDataType;

% For wav and flac files, specify the datatype based upon the 
% input data to preserve the incoming data and avoid
% conversions from inputDataType to channelDataType
filepath = multimedia.internal.io.FilePath(filename);
if any(strcmpi(filepath.Extension, {'wav','flac'}))
    switch inputDataType
        case 'uint8'
            dataType = 'int16'; % wav/flac device needs int16 minimum.
        case {'int16','int32','single','double'}
            dataType = inputDataType;
        otherwise
            dataType = channelDataType;
    end
end

end

function validateExtension (extension, filename)
% Validate that the filename has the correct extension.

if isempty(extension)
    error(message('MATLAB:audiovideo:audiowrite:noFileExtension', ...
        filename, sprintWriteableFileTypes));
end

import multimedia.internal.audio.file.PluginManager;
if ~any(strcmpi(extension, PluginManager.getInstance.WriteableFileTypes))
    error(message('MATLAB:audiovideo:audiowrite:invalidFileExtension',...
        filename, sprintWriteableFileTypes));
end

end

function [filename, exists] = validateFilename( filename )

filepath = multimedia.internal.io.FilePath(filename);

if isempty(filepath.Absolute)
    error(message('MATLAB:audiovideo:audiowrite:folderNotFound', ...
        filepath.ParentPath));
end
% Validate the file extension
validateExtension(filepath.Extension,filename)

if ~filepath.Writeable
    error(message('MATLAB:audiovideo:audiowrite:fileNotWritable', filename));
end

filename = filepath.Absolute;
exists = filepath.Exists;
end

function fileTypes = sprintWriteableFileTypes
% Print the list of writeable file types to a string
import multimedia.internal.audio.file.PluginManager;
fileTypes = cellfun(@(x) sprintf('\t.%s\n',x),...
    PluginManager.getInstance.WriteableFileTypes, ...
    'UniformOutput', false);

% expand into a character array
fileTypes = [fileTypes{:}]; 

% remove trailing newline ('\n') character.
fileTypes = fileTypes(1:end-1);
end

function warnIfMp4Metadata(props)
% Warn when writing Mp4 metadata on Mac, which is currently unsupported
if ~ismac
    % Metadata s
    return
end

if isempty(props.Title) && isempty(props.Artist) && isempty(props.Comment)
    % If no metadata is being written, do not warn
    return
end

% Warn only for m4a and mp4 files
filepath = multimedia.internal.io.FilePath(props.filename);
if any(strcmpi(filepath.Extension, {'m4a','mp4'}))
    warning(message('MATLAB:audiovideo:audiowrite:metadataMp4UnsupportedMac'));
end
end

function props = validateProperty(propName, props,channelProp, channel, defaultValue)
% A property (propName) is only valid if its associated
% channel property (channelProp) is available on the channel.
if ~isempty(props.(propName)) && ~ismember(channelProp, properties(channel))

    [~,~,extension] = fileparts(props.filename);
    error(message(...
        'MATLAB:audiovideo:audiowrite:unsupportedParameter',...
        propName,extension));
end

% Property is valid.  Set the default value if it has not been
% specified by the user
if isempty(props.(propName))
    props.(propName) = defaultValue;
end

end

%function y = clipInputData(y, fileDataType)
%if ~ismember(class(y),{'double','single'})
%    return; % data does not need to be clipped
%end

%if ismember(fileDataType, {'single','double'})
    % do not clip when writing out 'single' or 'double' to the file
    % to preserve full fidelity of the signal
 %   return; 
%end

% Data is input as a 'single' or 'double', and being written to an integer
% format.  Clip the range [-1 +1].
%if any(y(:) < -1 | y(:) > +1)
%    warning(message('MATLAB:audiovideo:audiowrite:dataClipped')); 
%end
%y(y < -1) = -1;
%y(y > +1) = +1;

%end

function handleException(exception, props, fileExisted)
import multimedia.internal.audio.file.PluginManager;

newexception = PluginManager.convertPluginException(exception, ...
                                           'MATLAB:audiovideo:audiowrite');
    
filepath = multimedia.internal.io.FilePath(props.filename);

if ~fileExisted && filepath.Exists
    % delete new files on error
    delete(props.filename);
end

throwAsCaller(newexception);

end
