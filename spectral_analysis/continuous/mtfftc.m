function [J,freq]=mtfftc(data,tapers,nfft,Fs,usewavelet)
% Multi-taper fourier transform - continuous data
%
% Usage:
% J=mtfftc(data,tapers,nfft,Fs) - all arguments required
% Input: 
%       data (in form samples x channels/trials or a single vector) 
%       tapers (precalculated tapers from dpss) 
%       nfft (length of padded data)
%       Fs   (sampling frequency)
%                                   
% Output:
%       J (fft in form frequency index x taper index x channels/trials)
if nargin < 4; error('Need all input arguments'); end
if nargin < 5, usewavelet=false; end 

data = change_row_to_column(data);
[NC,C]=size(data); % size of data
[NK K]=size(tapers); % size of tapers
if NK~=NC; error('length of tapers is incompatible with length of data'); end;
tapers=tapers(:,:,ones(1,C)); % add channel indices to tapers
data=data(:,:,ones(1,K)); % add taper indices to data
data=permute(data,[1 3 2]); % reshape data to get dimensions to match those of tapers
data_proj=data.*tapers; % product of data with tapers
if usewavelet
    s0=2*(1/Fs);
    ds = 0.4875;
    sig=struct();
    sig.val = data_proj;
    sig.period = 1/Fs;
    NbSc=fix(log2(length(sig.val)*(1/Fs)/s0)/ds);
    scales = struct('s0',s0,'ds',ds,'nb',NbSc);
    for i = 1:size(sig.val,2)
        temp=cwtft(sig,'scales',scales);   % fft of projected data
        freq = temp.frequencies;
        J(:,i,:) = temp.cfs;
    end
else
    J=fft(data_proj,nfft)/Fs;   % fft of projected data
end
