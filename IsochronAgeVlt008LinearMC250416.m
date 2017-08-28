

% This File Compilates measurement data, assembles constants, and
% caluclate the isochron age by calling york.m (by Balco and Rovey, 2008)
% and billipse.m (by Balco and Rovey, 2008)
% Error caclulation done by Monte Carlo simulation ?


clear all; close all;

%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Input Data for Isochron Dating of Terrace VLT008 ref. to 012
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


% Measurements of 10Be, absolute error, 26Al, and absolute error

%yearly DK blank
                  data_XXX = [%128985	21669	919376	169523
                                1004223	17091	2396862	145176
                                177737	5389	673464	65106
                                91573	4808	577993	49902
                                108019	8931	512452	51601
                                168400	9912	712168	53685
                                1014001	 16619	2422200	226320
                                93076	4118	351921	43458];
                                %154758	5301	1102775	58588];
                            
                             
%Individual DK blank
%                   data_XXX = [128985	21669	737449	158990
%                                 1004223	17091	2355536	144560
%                                 177737	5389	632206	63724
%                                 91573	4808	537990	48197
%                                 108019	8931	471191	49847
%                                 168400	9912	670331	51953
%                                 1014001	16619	2382446	225954
%                                 93077	4118	353984	44042
%                                 154758	5301	1106028	59662];
                     
                            
                            
%yearly AV blank
%                   data_XXX = [128985	21669	994071	137117
%                                 1004223	17091	2413829	143399
%                                 177737	5389	690404	61055
%                                 91573	4808	594417	44831
%                                 108019	8931	529393	46385
%                                 168400	9912	729345	48545
%                                 1014001	16619	2438522	225269
%                                 93077	4118	336283	43497
%                                 154758	5301	1078121	58661];
                     
                     

                
                       
% Sampling depth [cm] and terrace density [g/cm3]

    depth = 1210;
    deldepth = 100;
    
    rho = 2.00;
    delrho = 0.1;


% First guess of denudation rate (cm/years)

    eros_XXX = [0.0015
                %0.0015
                %0.0015
                0.0015
                0.0015
                0.0015
                0.0015
                0.0015
                0.0015];
    
    eros.ue = 0;
    
    
   % First guess of Beat value as in Granger adn Muzikar, 2001
       beta_XXX = [1
           1
            1
            1
            1
            1
            1
            1
           1];
    
% Be-10 measurements and uncertainties (atoms/g)
    data.x = data_XXX(:,1);
    data.dx = data_XXX(:,2);
% Al-26 measurements and uncertainties (atoms/g)
    data.y = data_XXX(:,3);
    data.dy = data_XXX(:,4);
    
% Denudation rates
    eros.z = eros_XXX(:,1);

% Beta value
    beta.z = beta_XXX(:,1);
      
 % 10Be and 26Al Production rate at sample surfac      
    data.P100sp = 5.275; 
    data.dP100sp = 0.222;
    data.P100sm = 0.028;
    data.P100fm = 0.037;
    
    data.P260sp = 35.444; 
    data.dP260sp = 1.489;
    data.P260sm = 0.343;
    data.P260fm = 0.357;

       
 % 10Be amd 26Al Production rate in catchement area
    data.P10csp = 7.079; 
    data.dP10csp = 0.297
    data.P10csm = 0.032;
    data.P10cfm = 0.038;
    
    data.P26csp = 47.407; 
    data.dP26csp = 1.991;
    data.P26csm = 0.396;
    data.P26cfm = 0.363;
    
    
    
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% General Input Data for Isochron Dating
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   
  % Decay constants
    l10 = 4.99746e-07;                 
    dell10 = 0.43e-08;
    l26 = 9.8300e-07;
    dell26 = 2.5000e-08;
    
   
  % Absorption length in g/cm2
    absorp(1) = 157;                  
    absorp(2) = 1030;
    absorp(3) = 160;
    absorp(4) = 3000;
    absorp(5) = 100;
    absorp(6) = 1520;
    absorp(7) = 7600;
    

  % 10Be Production parameters for sample location
    a(1) = 1 * data.P100sp ;            
    a(2) = 0.845 * data.P100sm;
    a(3) = -0.05 * data.P100sm;
    a(4) = 0.205 * data.P100sm;
    a(5) = 0.01 * data.P100fm;
    a(6) = 0.615 * data.P100fm;
    a(7) = 0.375 * data.P100fm;
    
  % 26Al Production parameters for sample location
    b(1) = 1 * data.P260sp ;           
    b(2) = 0.845 * data.P260sm;
    b(3) = -0.05 * data.P260sm;
    b(4) = 0.205 * data.P260sm;
    b(5) = 0.01 * data.P260fm;
    b(6) = 0.615 * data.P260fm;
    b(7) = 0.375 * data.P260fm;
    
  % 10Be Production parameters for sample depth
    e(1) = 1 * data.P100sp  * exp(- rho * depth / (absorp(1)));          
    e(2) = 0.845 * data.P100sm * exp(- rho * depth / (absorp(2)));
    e(3) = -0.05 * data.P100sm * exp(- rho * depth / (absorp(3)));
    e(4) = 0.205 * data.P100sm * exp(- rho * depth / (absorp(4)));
    e(5) = 0.01 * data.P100fm * exp(- rho * depth / (absorp(5)));
    e(6) = 0.615 * data.P100fm * exp(- rho * depth / (absorp(6)));
    e(7) = 0.375 * data.P100fm * exp(- rho * depth / (absorp(7)));
    
  % 26Al Production parameters for sample depth
    f(1) = 1 * data.P260sp  * exp(- rho * depth / (absorp(1)));           
    f(2) = 0.845 * data.P260sm * exp(- rho * depth / (absorp(2)));
    f(3) = -0.05 * data.P260sm * exp(- rho * depth / (absorp(3)));
    f(4) = 0.205 * data.P260sm * exp(- rho * depth / (absorp(4)));
    f(5) = 0.01 * data.P260fm * exp(- rho * depth / (absorp(5)));
    f(6) = 0.615 * data.P260fm * exp(- rho * depth / (absorp(6)));
    f(7) = 0.375 * data.P260fm * exp(- rho * depth / (absorp(7)));
    
    
  % 10Be Production parameters for catchment area   
    c(1) = 1 * data.P10csp ;            
    c(2) = 0.845 * data.P10csm;
    c(3) = -0.05 * data.P10csm;
    c(4) = 0.205 * data.P10csm;
    c(5) = 0.01 * data.P10cfm;
    c(6) = 0.615 * data.P10cfm;
    c(7) = 0.375 * data.P10cfm;

   % 26Al Production parameters for catchment area
    d(1) = 1 * data.P26csp ;           
    d(2) = 0.845 * data.P26csm;
    d(3) = -0.05 * data.P26csm;
    d(4) = 0.205 * data.P26csm;
    d(5) = 0.01 * data.P26cfm;
    d(6) = 0.615 * data.P26cfm;
    d(7) = 0.375 * data.P26cfm;
    

    % Production ratio at sample location
    data.Rp = (data.P260sp + data.P260sm + data.P260fm)/(data.P100sp + data.P100sm + data.P100fm); 
    
    % Production ratio at sample depth
    
    P10dRp = 0;
    P26dRp = 0;
    
    for r = 1:7  
        P10dRp = P10dRp + e(r);
        P26dRp = P26dRp + f(r);
    end  
    
    data.dRp = P26dRp / P10dRp;
    
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 % Start Calculation of Age
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    sdata.x = data.x;
    sdata.dx = data.dx;
    sdata.y = data.y;
    sdata.dy = data.dy;


% Take an initial guess at the age. 
 
    [initslope, initslopedel, intercept, interceptdel] = york(data)
    initage = (-log(initslope./data.Rp))./(l26-l10);

    disp(['Initial slope guess is ' num2str(initslope) ]);
    disp(['Initial slope error is ' num2str(initslopedel) ]);
    disp(['Initial age guess is ' num2str(initage./1e6) ' Myr']);
    
    
       
%  Main iteration loop

    tol = 10; % Tolerance in years. This is like 1e-3 accuracy. 


    oldage = initage;
    oldintercept = intercept;
    oldinitslope = initslope;
    outcount = 0;
            

    
    
 while 1    % loop: deals with the postdepositional production
               % Figure the postdepositional production
        
        Npb26.y = 0;
        Npb10.x = 0;     


        g = data.dRp * ((l10 * (1 - (exp(- l26 * oldage)))) ) ./ ((l26 * (1 - (exp(- l10 * oldage)))));
        
        Npb10.x = oldintercept / (g - oldinitslope);
        Npb26.y = Npb10.x * g;
        %Npb26.y = Npb10.x * data.Rp;
                 
        Neros10.x = (data.x - Npb10.x) / (exp(- l10 * oldage));
        Neros26.y = (data.y - Npb26.y) / (exp(- l26 * oldage));
        
        
         % Iterative solution of cosmogenic denudation rate for each
         % sample`
         
            Neros10calc.x = 0; 
        
         for s = 1:7  
                    Neros10calc.x = Neros10calc.x + (c(s) ./ (l10 + (eros.z * 2.7 / (absorp(s))))); 
         end
        
        
           while abs ((Neros10.x ./ (Neros10calc.x)) - 1) > 0.0001
              eros.z = eros.z - eros.z .* ((Neros10.x - Neros10calc.x ) ./ Neros10.x);
                               
               
              Neros10calc.x = 0;

                for t = 1:7  
                    Neros10calc.x = Neros10calc.x + (c(t) ./ (l10 + (eros.z * 2.7 / (absorp(t))))); 
                end
                       
           end  
           
           
         % Determination of linearization factor?????????
           
               h.z = 0;
               P10eros.z = 0;
               P26eros.z = 0;
               
         for u = 1:7
               P10eros.z = P10eros.z + c(u) ./ (l10 + ((eros.z * 2.7) / absorp(u)));
               P26eros.z = P26eros.z + d(u) ./ (l26 + ((eros.z * 2.7) / absorp(u)));  
         end
         
         Rinh.z = P26eros.z ./ P10eros.z; 
         h.z = Rinh.z ./ data.Rp;
         
        corrdata.x = h.z .* data.x;
        corrdata.y = data.y;
        
        sdata.x = corrdata.x;
        sdata.y = corrdata.y;
        
        resetN10 =  (data.x - Npb10.x) / (exp(- l10 * initage));
        resetN26 = (data.y - Npb26.y) / (exp(- l26 * initage));
		
      
        
        % Figure the new age;
        [newR,newdelR,newa,newdela,newdiag] = york(sdata);   
        newage = (-log(newR/data.Rp))./(l26-l10);
        

    
     
    % End the loop if adequately converged
    outcount = outcount + 1;
    
    if abs(newage - oldage) < tol;
       disp(['loop done -- ' int2str(outcount) ' iterations']);
         break;
       end;
    
     
    if outcount > 20;
       disp(['loop bailed after 20 iterations']);
         break;
       end;

    oldage = newage;
    oldinitslope = newR;
    oldintercept = newa;
  
    beta.z = Npb10.x ./  (Neros10calc.x * (exp(- l10 * newage)));
    %beta.z = Npb10.x ./  (data.x - Npb10.x);
    
 end;
 
 
     disp(['Final age guess is ' num2str(newage./1e6) ' Myr']);
     disp(['Final slope is ' num2str(newR) ]);
     disp(['Final slope error is ' num2str(newdelR) ]);

     %disp(['beta is' num2str(beta) 'mm/yr']);
     %disp(['denudation rate is ' num2str(eros * 10) 'mm/yr']);

													
    
  
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 % Figure 1. The isochron diagram.
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    figure(1); clf; 

    % Data as measured -- blue ellipses
    plot(data.x,data.y,'b.');hold on;
    for a = 1:length(data.x);
        billipse(data.x(a),data.dx(a),data.y(a),data.dy(a),1,'b');
    end;
        
    % Data corrected for postdepostional production -- black ellipses
    plot(corrdata.x,corrdata.y,'k.');
    for a = 1:length(corrdata.x);
            billipse(corrdata.x(a),data.dx(a),corrdata.y(a),data.dy(a),1,'k');
    end;
    
    
        % Data corrected for postdepositional production and un-decayed
    plot(resetN10,resetN26,'r.');
    for a = 1:length(resetN10);       
           billipse(resetN10(a),data.dx(a),resetN26(a),data.dy(a),1,'r');
    end;

    
    % Line at slope Rp
    xx = [0 1.5e6]; yy = xx.*data.Rp;
    plot(xx,yy,'k');
    
    % Isochron
    xx = [0 1.5e6];
    yy = newa + xx.*newR;
    plot(xx,yy,'g');
    
    % Text
    %text(-5.5e5,2.5e5,agestr);
    xlabel('[Be-10] (atoms/g)');
    ylabel('[Al-26] (atoms/g)');
    title('Isochron diagram');
 
    % clean up
    axis([0 15e5 0 10e6]);
    
    % End of Figure 1
% 
% %xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


    
%    
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
% Start Monte Carlo Calculation of Age and Denudation Rate Error
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    
    % ********
    % *** START USER DEFINED PARAMETES
    % ********
        % number of monte carlo simulations to try
            iter = 50000;
		
        % Define the minimum and maximum value for each free parameter
        
        newR_min = newR - newdelR;
        newR_max = newR + newdelR;
        
        data.x_min = data.x - data.dx;
        data.x_max = data.x + data.dx;
        
        data.y_min = data.y - data.dy;
        data.y_max = data.y + data.dy;
        
        data.P100sp_min = data.P100sp - data.dP100sp;   %[atms/gr/yr]  minumum production rate at sample locatioin for nucleons
        data.P100sp_max = data.P100sp + data.dP100sp;   %[atms/gr/yr]  maximum production rate at sample locatioin for nucleons
        
        data.P260sp_min = data.P260sp - data.dP260sp;   %[atms/gr/yr]  minumum production rate at sample locatioin for nucleons
        data.P260sp_max = data.P260sp + data.dP260sp;   %[atms/gr/yr]  maximum production rate at sample locatioin for nucleons
	 
    	rho_min = rho - delrho;				% [gr/cm3] minimum density of soil profile
        rho_max = rho + delrho;     		% [gr/cm3] max density of soil profile
	
        depth_min = depth - deldepth;
        depth_max = depth + deldepth;
        
        l10_min = l10 - dell10;          	% [1/yr] minimum decay constant
        l10_max = l10 + dell10;          	% [1/yr] maximum decay constant 
        
        l26_min = l26 - dell26;          	% [1/yr] minimum decay constant
        l26_max = l26 + dell26;          	% [1/yr] maximum decay constant 
  
        
        sum_newRMC = 0;
        sum_del_newRMC = 0;
        sum_newageMC = 0;
		sum_del_newageMC = 0;
        sum_erosMC.z = 0;
        sum_del_erosMC.z = 0;
        sum_Npb10MC.x = 0;
        sum_del_Npb10MC.x = 0;
        sum_Neros10MC.x = 0;
        sum_del_Neros10MC.x = 0;
                  
    % ********
    % *** END USER DEFINED PARAMETES
    % ********	
    
        % *** Calculate some initial values needed in monte carlo loop
        newRMC_output = zeros(iter);
        newageMC_output = zeros(iter);
        erosMC_output = zeros(iter);
	
    % ********
    % *** START MONTE CARLO SIMULATION
    % ********	
    
    
for k = 1:iter
        %display('Iteration number:  ') 
        k;
        
        
        % Denudation rates
        initerosMC_XXX = [0.0015
                0.0015
                0.0015
                0.0015
                0.0015
                0.0015
                0.0015];
                    
        initerosMC.ue = 0;
        
        initerosMC.z = initerosMC_XXX(:,1);

        
	% *** Calculate random values between the min and max value for each free parameters for this iteration
		
        initslopeMC = rand * (newR_max - newR_min) + newR_min;
        data.x = rand * (data.x_max - data.x_min) + data.x_min;
        data.y = rand * (data.y_max - data.y_min) + data.y_min;
        data.P100sp = rand * (data.P100sp_max - data.P100sp_min) + data.P100sp_min;
		data.P260sp = rand * (data.P260sp_max - data.P260sp_min) + data.P260sp_min;
		rho = rand * (rho_max - rho_min) + rho_min;
        depth = rand * (depth_max - depth_min) + depth_min;
		l10 = rand * (l10_max - l10_min) + l10_min;
		l26 = rand * (l26_max - l26_min) + l26_min;
		
		% *** Done with assigning radom values
        
        % Check if ther are negative values.
        
         if data.x <= 0;
             break;
         end
                 
         if data.y <= 0;
             break;
         end
        
         
        
	
		% Take an initial guess at the age. 
 
            [initslopeMC, initslopdelMC, interceptMC, interceptdelMC] = york(data);
            initageMC = (-log(initslopeMC./data.Rp))./(l26-l10);

            %disp(['Age guess is ' num2str(age./1e6) ' Myr']);
    
        
            
                
%  Main iteration loop
 
    tol = 10; % Tolerance in years. This is like 1e-3 accuracy. 
 
 
    oldageMC = initageMC;
    erosMC.z = initerosMC.z;
    oldinterceptMC = interceptMC;
    oldinitslopeMC = initslopeMC;
    outcount = 0;
            


    
 while 1    % loop: deals with the postdepositional production
               % Figure the postdepositional production
        
        Npb26MC.y = 0;
        Npb10MC.x = 0;     
   
            

        g = data.dRp * ((l10 * (1 - (exp(- l26 * oldageMC)))) ) ./ ((l26 * (1 - (exp(- l10 * oldageMC)))));
        
        Npb10MC.x = oldinterceptMC / (g - oldinitslopeMC);
        Npb26MC.y = Npb10MC.x * g;
                 
        Neros10MC.x = (data.x - Npb10MC.x) / (exp(- l10 * oldageMC));
        Neros26MC.y = (data.y - Npb26MC.y) / (exp(- l26 * oldageMC));
        

         % Iterative solution of cosmogenic denudation rate for each
         % sample`
         
            Neros10calcMC.x = 0; 
        
         for s = 1:7  
                    Neros10calcMC.x = Neros10calcMC.x + (c(s) ./ (l10 + (erosMC.z * 2.7 / (absorp(s))))); 
         end
        
        
           while abs ((Neros10MC.x ./ (Neros10calcMC.x)) - 1) > 0.0001
              erosMC.z = erosMC.z - erosMC.z .* ((Neros10MC.x - Neros10calcMC.x ) ./ Neros10MC.x);
    
                Neros10calcMC.x = 0;
 
                for t = 1:7  
                    Neros10calcMC.x = Neros10calcMC.x + (c(t) ./ (l10 + (erosMC.z * 2.7 / (absorp(t))))); 
                end
        
           end  

           
       
    % Determination of linearization factor?????????
           
               h.z = 0;
               P10erosMC.z = 0;
               P26erosMC.z = 0;
               
         for u = 1:7
               P10erosMC.z = P10erosMC.z + c(u) ./ (l10 + ((erosMC.z * 2.7) / absorp(u)));
               P26erosMC.z = P26erosMC.z + d(u) ./ (l26 + ((erosMC.z * 2.7) / absorp(u)));  
         end
         
         Rinh.z = P26erosMC.z ./ P10erosMC.z; 
         h.z = Rinh.z ./ data.Rp;
         
        corrdata.x = h.z .* data.x;
        corrdata.y = data.y;
        
        sdata.x = corrdata.x;
        sdata.y = corrdata.y;
        
        resetN10 =  (data.x - Npb10MC.x) / (exp(- l10 * oldageMC));
        resetN26 = (data.y - Npb26MC.y) / (exp(- l26 * oldageMC));  
    
               
       if sdata.x <= 0;
             break;
       end
       
       
       % Figure the new age;
        [newRMC,newdelRMC,newaMC,newdelaMC,newdiagMC] = york(sdata);   
        newageMC = (-log(newRMC/data.Rp))./(l26-l10);
       
    
    % End the loop if adequately converged
    outcount = outcount + 1;

    
    if abs(newageMC - oldageMC) < tol;
       %disp(['loop done -- ' int2str(outcount) ' iterations']);
         break;
       end;
    
     
    if outcount > 20;
       %disp(['loop bailed after 20 iterations']);
         break;
       end;
 
    oldageMC = newageMC;
    oldinitslopeMC = newRMC;
    oldinterceptMC = newaMC;
    
   end;
 

 	% *** Assign calcualted output values to a output vectors for analysis after the monte carlo loop is done
        
        newRMC_output(k) = newRMC;
		newageMC_output(k) = newageMC;
		%eros_output(k).z = eros.z;
	
	    sum_newRMC = sum_newRMC + newRMC;
        sum_newageMC = sum_newageMC + newageMC;
        sum_erosMC.z = sum_erosMC.z + erosMC.z;
        sum_Npb10MC.x = sum_Npb10MC.x + Npb10MC.x;
        sum_Neros10MC.x = sum_Neros10MC.x + Neros10MC.x;
	
end;
  
     
  %calculate meag age and mean denudation rate
  
      mean_newRMC = sum_newRMC / k;  
      mean_newageMC = sum_newageMC / k;
      mean_erosMC.z = sum_erosMC.z / k  ;
      mean_Npb10MC.x = sum_Npb10MC.x / k;
      mean_Neros10MC.x = sum_Neros10MC.x / k;
      
      for j = 1:iter
        sum_del_newRMC = sum_del_newRMC + (initslopeMC - mean_newRMC).^ 2;
        sum_del_newageMC = sum_del_newageMC + (initageMC - mean_newageMC) .^ 2;
        sum_del_erosMC.z = sum_del_erosMC.z + (initerosMC.z - mean_erosMC.z) .^ 2;
        sum_del_Npb10MC.x = sum_del_Npb10MC.x + (Npb10.x - mean_Npb10MC.x) .^ 2;
        sum_del_Neros10MC.x = sum_del_Neros10MC.x + (Neros10.x - mean_Neros10MC.x) .^ 2;
      end  
      
      del_newRMC = sqrt ((sum_del_newRMC) / (k - 1));  
      del_newageMC = sqrt ((sum_del_newageMC) / (k - 1));
      del_erosMC.z = sqrt ((sum_del_erosMC.z) / (k - 1));
      del_Npb10MC.x = sqrt ((sum_del_Npb10MC.x) / (k - 1));
      del_Neros10MC.x = sqrt ((sum_del_Neros10MC.x) / (k - 1));
      
      final_del_newRMC = del_newRMC ./ mean_newRMC * 100;
      final_del_newageMC = del_newageMC ./ mean_newageMC * 100;
      final_del_erosMC.z = del_erosMC.z ./ mean_erosMC.z * 100;
      final_del_Npb10MC.x  = del_Npb10MC.x  ./ mean_Npb10MC.x  * 100;
      final_del_Neros10MC.x = del_Neros10MC.x./ mean_Neros10MC.x* 100;
      
      den_rate.z  =  eros.z .* 10;
      del_den_rate.z  = eros.z .* final_del_erosMC.z / 100 * 10;
    
      del_Npb10MC.x  = Npb10.x  .* final_del_Npb10MC.x / 100;
      del_Neros10MC.x   = Neros10.x  .* final_del_Neros10MC.x / 100;

display('**********************************') 
display('Monte Carlo Iterations Completed  ') 
display('**********************************') 
	
% ********
% *** CALCULATE VALUES AND ERRORS OF AGE AND DENUDATION RATE
% ********		
disp(['Slope is ' num2str(newR) ]);
disp(['Slope error is ' num2str(newR * final_del_newRMC / 100) ]);

disp(['Age is ' num2str(newage) 'yr']);
disp(['Error is ' num2str(newage * final_del_newageMC / 100) 'yr']);

disp(['Beta is']);
beta.z

disp(['Denudation rate is mm/yr']);
den_rate.z
disp(['Error is mm/yr']);
del_den_rate.z

disp(['Post depos. concentration ']);
Npb10.x
disp(['Error ']);
del_Npb10MC.x

disp(['inherited concentration']);
Neros10.x
disp(['Error ']);
del_Neros10MC.x


display('**********************************') 
	
	
    
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx



