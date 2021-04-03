%% To convert Time series text file to matlab Readable format
clear;
fid=fopen('1752318_TS5_1_2.txt');    % Load a file
% delete outfile.txt
nr=0;                       % Number of records
samp_int=1/Fs;              % Fs sampling freqency
fidout=fopen('output_1_2.txt','w');   % output file
 
tic
while ~feof(fid) 
    A=[];
    block_data=[];
    Tint=[];
    nr=nr+1;
    block_data=textscan(fid, repmat('%d',[1,6]), 'HeaderLines',5,'CollectOutput',1); % To remove header files from the text file
    A=cell2mat(block_data);
    Tint=samp_int+nr : samp_int : nr+1; Tint=Tint';
    nr
    if(length(Tint)==length(A(:,1)))
        for i=1:length(A(:,1))
            fprintf(fidout,'%12.4f %8d %8d %8d %8d %8d %8d\n',Tint(i),A(i,:));
        end
    end
end
disp(['Total # records processed :',num2str(nr),' sec'])
fclose(fid);
fclose(fidout);
%% To plot time series
A=load('output_56.txt');
[row,col]=size(A);
for i=2:col
    figure(i)
    plot(A(:,1),A(:,i),'.-');
    title(['Channel: ',num2str(i-1)]);
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
    fig_fname=['Channel',num2str(i-1),'.jpg'];
    saveas(gcf,fig_fname);
end

%% power spectral density plot of six channels
fs=15;
A=load('output_1_2.txt');
name = 'xyz';
for col_idx = 2:4
    f = figure;
    for offset = 0:3:3
        data=A(:,col_idx+offset);  %column data
        nfft=1024;                 % nonequispaced fast fourier transform in multiple of 2
        Y=fft(data,nfft);          %fourier transform that carries complex values
        Y=Y(1:nfft/2+1);
        f=((0:1/nfft:1/2)*fs).';    % Frequency bins
        magnitude=(1/(fs*nfft))*(abs(Y)).^2;
        magnitude(2:end-1)= 2*magnitude(2:end-1);
        dB_mag = magnitude;
        %dB_mag=10*log10(magnitude);  % to convert magnitude in log
        if offset == 0
            subplot(1,2,1);
        else
            subplot(1,2,2);
        end
        loglog(f,dB_mag);
        pbaspect([1 0.75 1]);
        grid on;
        set(gcf, 'Position', [50, 50, 1000, 600]);
        title('Magnitude response of signal');
        if offset == 0
            ylabel(['Magnitude of E_',name(col_idx-1),'(V/sqrt(Hz))']);
        else
            ylabel(['Magnitude of H_',name(col_idx-1),'(nT/sqrt(Hz)']);
        end
        xlabel('Frequency(Hz)')
        hold on
    end
    hold off
end