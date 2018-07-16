classdef PL
% PL - Phase Library
%
% Author: Ryan Y
% Library of phase measures one can unleash that rival or are better than
% standard coherence.
%
% TODO VERIFY WPLI AND DUAL-F Coherence

  methods(Static) % Given this class is a static library, just has this methods section

    function [plv,pli,wpli,imc,du_imc,du_wpli] ...
        = computemeasures(J1,J2,C12,varargin)
    %PHASEAMPDOUBLE_COHERENCE Returns the phase, amplitude, and cross
    %coherency given two fourier series.
    
      do_dual = true;
      phase = false;
      optlistassign(who,varargin{:});
      sqmean=@(x,y) squeeze(mean(x,y));
      
      if phase == true
        J1=J1./abs(J1);
        J2=J2./abs(J2);
        J12 = J1.*J2;
        C12 = sqmean(J12,2)/sqrt(sqmean(J1,1).*sqmean(J2,2));
      else
        J12 = J1.*conj(J2);
      end
      
      %% Non-dual measures
      
      plv = PL.plv(J12);
      pli = PL.pli(J12);
      wpli = PL.wpli(J1,J2,J12);
      imc = PL.imc(C12);
      
      %% Dual measures
      if do_dual && nargout > 4
        du_imc=[];
        du_wpli=[];
        if do_dual
          du_imc       = PL.du_imc(J1,J2);      % Dual imaginary coherence
          %du_wpli      = PL.du_wpli(J1,J2);     % Dual weighted phase lag index
        end
      end

    end
    % --------------------------------------------------------------------  
    % ---------------------Single Coherent Measures-----------------------  
    % --------------------------------------------------------------------  
    function plv = plv(J12,dim)
    % PLV computes the phase-locking value
    % Equivalent calc abs(mean( S12/abs(S12) ,2))
      if nargin == 1 ; dim = 2; end
      plv = abs(mean( J12./abs(J12),dim));
      % Alternative way to calculate the same number
      %plv = sqrt(mean( cos(angle(J12)), dim).^2 + mean(sin(angle(J12)),dim).^2); % compare this to the above!
      %plv = abs(mean( exp(angle(J12)*1j) , dim ));
    end
    % --------------------------------------------------------------------
    function pli = pli(J12,dim)
    % CALC_SPI computes the phase-locking value
      if nargin==1; dim = 2; end
      pli = abs(mean( sign(imag(J12)) ,dim));
    end
    % --------------------------------------------------------------------
    function wpli = wpli(J1,J2,J12,varargin) %#ok<INUSL>
    % CALC_WSPI computes the weighted phase locking index
    %
    % Version 2 - replacing with my own
    % 
    % I tried to replicate the equation in the paper for the debiased but I
    % might be thinking about X^(1) and X^(2) wrongly. I'll need to read
    % much deeper to be sure I've got it.
    
      dim = 2;
      func=@sum;
      debias = true;
      optlistassign(who,varargin{:});
      
      % The method used here to compute the expected value of E{ imag(Xp) imag(Xq) }
      % over tapers or trials who are each samples of some true underlying X cross-spectra
      % is essentially taking the different of squares to find the squared wpli...the squared
      % wpli is an debiased estimator
      
      J = imag(J12);
      a = abs(func(J,dim));
      J(isnan(J))=0;
      b = func(abs(J),dim);
      if debias
       c = func(J.^2,dim);
       wpli = (a.^2 - c) ./ (b.^2 - c);
      else
       wpli = a./b;
      end

      %if debias
      %  improd = imag(J1).*imag(J2);
      %  wpli = mean(improd,2)./mean(abs(improd),2);
      %else
      %  wpli = abs(mean(imag(J12),2))./mean(abs(imag(J12)),2);
      %end
      
    end
    % --------------------------------------------------------------------  
    function ImC = imc(C12)
      % CALC_IMC computes imaginary coherence
        ImC = abs(imag(C12));
    end
    % --------------------------------------------------------------------  
    % -----------------------Dual Coherent Measures-----------------------  
    % --------------------------------------------------------------------  
    function [Full,Imaginary] = du_C(J1,J2,varargin)
      
      % Process the parameters
      dim = 2;
      gaussfilt = 0; %3
      S1 = []; S2 = [];
      type='normal';
      optlistassign(who,varargin{:});
      
      % If user did not provide power spectra 
      if isempty(S1) || isempty(S2)
        S1=squeeze(mean(conj(J1).*J1,2));
        S2=squeeze(mean(conj(J2).*J2,2));
      end
      % Here, you can select phase-phase or phase--phase-amp coherence. 
      switch type
        case 'phase'
          J1 = J1./abs(J1);
          J2 = J2./abs(J2);
        case 'phase-phaseamp'
          J1 = J1./abs(J1);
      end
      
      % Repmat the spectral power densities over the entire extent of the
      % bispectral matrix
      Sx = repmat(S1,1,size(J1,1));
      Sy = repmat(S2,1,size(J1,1));
      
      % Number of trials or tapers
      K = size(J1,2);
      
      % Option to return full or imaginary component
      Jcross_im = imag(J1*conj(J2)')/K ;    % Mean of the imagianry cross
      Jcross_real = real(J1*conj(J2)')/K ;  % Mean of the real cross
      
      % Now get the coherence, imaginary coherence, and dual imaginary
      % coherence
      den=Sx.*Sy';
      Imaginary = (Jcross_im.^2)./den;
      Full = Imaginary + (Jcross_real.^2)./den;
      Imaginary=sqrt(Imaginary);
      Full=sqrt(Full);

      % If gaussfilt is given, then gaussian smooth
      if gaussfilt>0
        Imaginary = imgaussfilt(Imaginary,gaussfilt); 
        Full      = imgaussfilt(Full,gaussfilt); 
      end
      
    end
    % --------------------------------------------------------------------
    function du_wpli = du_wpli(J1,J2,varargin)
    % This measure doesn't exist in the lit, but I've generalized the
    % weighted phase lag index to the realm of dual coherence;
    
      % Process the parameters
      dim = 2;
      gaussfilt = 3;
      S1 = []; S2 = [];
      type='normal';
      optlistassign(who,varargin{:});
      
      % If user did not provide power spectra 
      if isempty(S1) || isempty(S2)
        S1=squeeze(mean(conj(J1).*J1,2));
        S2=squeeze(mean(conj(J2).*J2,2));
      end
      % Here, you can select phase-phase or phase--phase-amp coherence. 
      switch type
        case 'phase'
          J1 = J1./abs(J1);
          J2 = J2./abs(J2);
        case 'phase-phaseamp'
          J1 = J1./abs(J1);
      end
      
      % Number of trials or tapers
      K = size(J1,2);
      
      % Iterate through and mean the wpli measure
      num = zeros(size(J1,1),size(J2,1));
      den = zeros(size(J1,1),size(J2,1));
      for k = 1:size(J1,2)
        num = num + imag(J1(:,k) ) * imag(J2(:,k))'; 
        den = den + abs(imag(J1(:,k)) * imag(J2(:,k))');
      end
      num = num/K; % Complete the expectation function
      den = den/K; % Complete the expectation function
      du_wpli = num./den; % Create the 2D weighted phase lag index
      du_wpli(isnan(du_wpli))=0; % Replace any nans created by the denominator being zero
            
      % So now we have to sqrt
      %	du_wpli = sqrt(abs(du_wpli));
    end
    % --------------------------------------------------------------------
    function [out] = phasephase(J1, J2, n, m)
    %PHASEPHASE Description
    %	[OUT] = PHASEPHASE(J1, J2, n, m) implements measurement of n:m
    %	phase-phase coupling
      J1 = permute(J1, [1 4 2 3]);
      J2 = permute(J2, [4 1 2 3]);
      out = bsxfun(@minus, J1, J2);
      out = out./abs(out);
      out = squeeze(mean(out,3));
    end
    % --------------------------------------------------------------------
    function [inStruct] = battery_of_nm_coupling(inStruct, J1, J2)
    %BATTERY_OF_NM_COUPLING Runs a sequence of nm phase coupling measures
    %	[OUTSTRUCT] = BATTERY_OF_NM_COUPLING(INSTRUCT, J1, J2) 
    %	Runs a sequence of nm phase coupling measures and adds each of them to the
    %	inStruct.
        NM = [1,2; 1,3; 1, 5 ];
        for nm = NM'
            inStruct.(sprintf('nm%d%dpc',nm(1),nm(2))) = PL.phasephase(J1,J2,nm(1),nm(2));
        end
  end
%  -------------------------------------------------------------------

end
end

