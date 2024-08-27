h=[0.4 0 0.26 0 0 0.4 0 0.6 0 0.5];

waterfillings(200,h,16,2,1);
function output_power = waterfillings (Ptotal,noise,nSub,snr,noiseVar)
	% input is row
	Pnoise = [Ptotal ; -((snr * noiseVar)./abs(fft(noise,nSub)).^2)'];
	mask = logical(ones(1,length(Pnoise) -1)'); % to mask p only
	cPower = [1 repelem(-1 ,length(mask) -1 )]'; % just an intialization
	while ( sum(cPower < 0) ~= 0)
		cPower = calcWaterFilling(Pnoise(logical([1 ;mask])));
		mask(mask) = cPower(1:end -1) >= 0;
	end
	output_power = zeros(1,length(Pnoise) -1);
	output_power(mask') = cPower(1:end-1)';
	bar(-Pnoise(2:end)' + output_power) ; hold on ; bar(-Pnoise(2:end)');
	legend('Normalized noise power','Distributed Power');
	xlabel('Subchannels');ylabel('Power');
end
function output = calcWaterFilling (v)
	equation_matrix = [ones(1,length(v)-1) 0;
	eye(length(v)-1) repelem(-1,length(v) -1)'];
	output = equation_matrix \ v;
end