clear;
close all;
h=[0.4 0 0.26 0 0 0.4 0 0.6 0 0.5];
a=[2 1 4 0.5 1 2 1];
input_signal = randi([0 1] , 1 , 2*32*1000);
%input_signal = repmat([1 0 0 1] , 1,128);
%waterfilling(200,h,16,2,1);
%% work
main(input_signal,@(a,b,c) channels_awgn(a,b,c),...
1 , -20,15 , 'AWGN channel');
main(input_signal,@(a,b,c) channels_selective(a,b,c),...
h , -5,15 , 'Selective channel');

function main(input_signal,channel,h,noise_start,noise_end,figureTitle)
	BER = [];
	BER = [BER ; ofdm(input_signal,@(a) mapper_qpsk(a), 1,...
	8,16,...
	1,channel ,h,noise_start,noise_end,...
	[0].',[-1 1].',[1].')];
	BER = [BER ; ofdm(input_signal,@(a) mapper_qpsk(a), 3,...
	8,16,...
	1,channel ,h,noise_start,noise_end,...
	[0].',[-1 1].',[1].')];

	x = noise_start:noise_end;
	figure;
	semilogy(x,BER(1,:)); hold on; semilogy(x,BER(2,:));
	xlabel('Eb/No (dB)'); ylabel('BER') ; title(figureTitle);
	legend('QPSK uncoded' , 'QPSK coded');

end
function BER = ofdm(input_signal,mapper,elementSamples,...
	m,n,... %interleaver
	Eb,channel,h,noise_start,noise_end,...

	limits,levels,values)

	N0 = (Eb)./(10.^((noise_start:noise_end)./10));
	BER = zeros(1,length(N0));

	for i=1:length(N0)
		encoded_data = repeate_encoder(input_signal,elementSamples);
		interleaved_data = interleave(encoded_data,m*n,m,n);
		interleaved_data=encoded_data;
		mapped_data=mapper(interleaved_data);

		IFFT_signal = IFFT(mapped_data,64); % ifft with samples 64 in one row
		CEX_signal = add_cycilic_ext(IFFT_signal,16);

		channel_output = channel(CEX_signal,N0(i),h);
		% channel_output=CEX_signal;
		RX_data = reciever_ofdm(channel_output,16,limits,levels,values,h);
		deinterleaved_data = deinterleave( reshape(RX_data.',1,[]),m*n,m,n);
		%deinterleaved_data= reshape(RX_data.',1,[]);
		decoded_data = repeate_decoder(deinterleaved_data,elementSamples);
		BER(i) = sum(xor(decoded_data,input_signal))/length(input_signal);
	end

end

%% FFT , IFFT
function outputSignal = IFFT(input_signal,elementsCount)
	%input data is row
	% output is column
	%shaping is the length of data to be ifft together
	outputSignal = reshape(input_signal,elementsCount,[]).';
	%assumed length of signal is bigger than shaping
	for i=1:size(outputSignal,1)
		outputSignal(i,:) = ifft(outputSignal(i,:));
	end
end
function outputSignal = FFT(input_signal)
	% input is column (output of reciever)
	% output is column
	for i=1:rowsLength(input_signal)
		input_signal(i,:) = fft(input_signal(i,:));
	end
	outputSignal = reshape(input_signal.',1,[]);
end
%% CEX
function outputSignal = add_cycilic_ext(input_signal,n)
	outputSignal = [ input_signal(:,end-n+1:end) input_signal ];
end
function outputSignal = remove_cycilic_ext(input_signal,n)
	outputSignal = input_signal(:,n +1:end);
end
%% channel encoders
function [signal] = repeate_encoder(input_signal,count)
	% repeate every element in input_signal by count
	signal = repelem(input_signal,count);
end
function outputSignal = repeate_decoder(input_signal,count)
	%input signal is row
	outputSignal = reshape(input_signal,count,[]).';
	outputSignal = round(sum(outputSignal,2)/count); % soft decoding
	outputSignal = reshape(outputSignal,1,[]);
end
%% interleavers
function [signal] = interleave(input_signal,elementsCount,m,n)

	%signal is row
	signal = reshape(input_signal,elementsCount,[]).';
	for i =1:size(signal,1)
		signal(i,:) = matintrlv(signal(i,:),m,n);
	end
	signal = toRow(signal); %return row
end
function [signal] = deinterleave(input_signal,shaping,m,n)
	% row
	signal = reshape(input_signal,shaping,[]).';
	for i =1:(length(input_signal)/length(signal))
		signal(i,:) = matdeintrlv(signal(i,:),m,n);
	end
	signal = reshape(signal.',1,[]);
end
%% mappers
function [outputSignal] = mapper_qpsk(inputSignal)
	% signal *Alpha + c = maped value
	outputSignal = reshape(inputSignal,2,[]).';
	outputSignal = outputSignal .* [2 2j] + [-1 -1j];
	%sum Real and Imaginary parts
	outputSignal = toRow(sum(outputSignal,2));
end
%% channels
function outputSignal = channels_awgn(input_signal,N0,~)
	%input is column

	noise = @() sqrt(N0/2) .* randn(size(input_signal));
	n=noise();
	outputSignal = input_signal +n+n.*1j;
end
function outputSignal = channels_selective(input_signal,N0,h)
	%input is column with 80 sample
	noise = @() sqrt(N0/2) .* randn(size(input_signal));
	temp_signal = input_signal + noise() + 1j .* noise();
	outputSignal =zeros(size(input_signal));
	for i=1:size(outputSignal,1)
		o=conv(temp_signal(i,:),h);
		outputSignal(i,:) = o(1:size(input_signal,2));
	end
	%output is column
	%outputSignal = conv(reshape(outputSignal.',1,[]),h);
	%outputSignal = reshape(outputSignal(1:end -length(h)+1
	),length(input_signal(1,:)),[]).';
end
%% reciever
function outputSignal = reciever_ofdm(input_signal,cycilic_N,limits,levels,values,h)
	input_signal = remove_cycilic_ext(input_signal,cycilic_N);
	% input function is complex column
	if length(h) >1 %selective cahnnel mode
		temp_signal = zeros(size(input_signal));

		for i=1:size(input_signal,1)
		temp_signal(i,:) = fft(input_signal(i,:),64)./fft(h,64);
		end
		%outputSignal=temp_signal;
		outputSignal = reshape(temp_signal.',1,[]);
	else
		outputSignal = FFT(input_signal);
	end
outputSignal = reciever(outputSignal,1,limits,levels,values);
end
function [outputSignal] = reciever(input_signal , R , limits, levels,values)
	% all inputs are columns
	% output is column
	% quantize to get nearest point to the value on real and imag. axis.
	input_signal = input_signal ./ R;
	[~,real_part] = quantiz(real(input_signal), limits, levels);
	[~,imag_part] = quantiz(imag(input_signal), limits, levels);
	% first bit is the sign and rest of bits depends on the value on the
	% axis
	outputSignal = [real_part>0 demap(abs(real_part),unique(abs(levels)),values) ...
	imag_part>0 demap(abs(imag_part),unique(abs(levels)),values)];
end
%% utils
function [output] = demap(input , keys,values)
	% detecting values on axis (sign is already detected)
	if length(keys) == 1
		output = [];
	else
		[~, n] = size(values);
		output = zeros(length(input) , n);
		for i=1:length(keys)
			output = output + (input == keys(i)).* values(i,:);
		end
	end
end

function out= toRow(input_signal)
	out= reshape(input_signal.',1,[]);
end
function out = toColumns(input_signal,count)
	out=reshape(input_signal,count,[]).';
end

function [n] = rowsLength(input)
	% this function detects the number of rows in matrix
	n = numel(input)/numel(input(1,:));
end