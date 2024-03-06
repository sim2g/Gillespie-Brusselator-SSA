% -------------------------------------------------------------------
% ChemicalReaction.m
% The class ChemicalReaction is defined with it's properties and methods.
% -------------------------------------------------------------------
classdef ChemicalReaction
    % Inputting constant properties given in the assignment. These are hidden
    % as they are used for initialisation of working variables.
    properties (Constant, Hidden)
        CDX = [0.01 0.001]; % X reaction constants for part 1a. and 1b.
        CDY = [1.5 1.5]; % Y reaction constants for part 1a. and 1b.
        x0 = [300 3000]; % Initial X populations
        y0 = [350 3500]; % Initial Y populations
        t_min = 0; % Initial starting time
        blues = (1/256).*[ [0 142 204]; [101 147 245];[0 0 128]; [15 82 186];  [0 128 255]; ...
            [115 194 251]; [70 130 180]; [0 128 129]; [149 200 216]; [79 151 163]]; % Shades of blue for 10 realisations of the Y Populations
        reds = (1/256).*[[210 31 60]; [94 25 20]; [178 34 34]; [205 92 92]; [234 60 83]; ...
            [224 17 95]; [147 58 22];[180 55 87]; [141 2 31]; [66 13 9]]; % Shades of red for 10 realisations o the X Populations
        t_det = (0:0.01:4); % Linearly spaced vector of 400 elements to plot the deterministic solution for X and Y between 0 and 4 time units
        S0 = [0.1 1 10]; % Sizes of system in question 2.
        seed = 179; % Arbitrary random number seed.
    end
    % These are protected properties as they do not need to be accessed
    % directly from main.m
    properties (Access = protected)
        X % Number of Y molecules at a given time
        Y % Number of X molecules at a given time
        t % Current iteration time
        count % Number of tau-sized time-steps that have been evaluated
        t_sample % Values are recorded at every t_sample
        n_sample % The number of samples that have been taken
        n_max % The maximum number of tau-sized time-steps to be allowed.
        S % Working variable for S, assigned to one of the 3 values in S0 depending on which question part the object is representing.
        c % Working array for the 4 values of c in question 2 depending on the value of S.
        cdx % Working value for cdX in question 1 depending on if part a. or part b.
        cdy % Working value for cdY in question 1 depending on if part a. or part b.
        ydata % Array containing number of Y molecules at each t_sample value
        xdata % Array containing number of X molecules at each t_sample value.
        tdata % Array containing all of the t + tau time steps.
        h % Array containing the h values of the propensities for each object.
        a % The propensity of each reaction
        random % Array containing 2 random numbers
        tau % The time for the next reaction
        mu % Logical value which decides the next reaction
        a0 % The total propensity
        mean_ % Array containing the cumulative mean for each increasing t_sample
        standarddeviation % Array containing the cumulative standard deviation for each increasing t_sample
        skewness % Array containing the cumulative skewness for each t_sample as it increases
        Y_det % The deterministic value of the Y population at each value of t_det
        X_det % The deterministic value of the X population at each value of t_det
        Y_pred % The deterministic value of the Y population at each t_sample
        X_pred % The deterministic value of the X population at each t_sample
        xdatanorm % xdata minus the mean of xdata
        ydatanorm % ydata minus the mean of ydata
        pxx_y % Periodogram power spectral density (PSD) estimate of Y
        pxx_x % Periodogram power spectral density (PSD) estimate of X
        fx % The frequency vector in cycles / time unit of X population
        fy % The frequency vector in cycles / time unit of Y population
    end
    % t_max and t_delta are accessed directly from main.m in part 2c so
    % they need to be public
    properties (Access = public)
        t_max % The maximum time for the simulation
        t_delta % The sampling interval
    end
    % Defining public methods for the algorithm to be called in main.m.
    methods (Access = public)
        % The constructor function is defined which creates an object for
        % each realisation or group of realisations. m is the question
        % number as given in the assignment and n is the index of cdX,
        % cdY, x0 or y0, or the index of S0, selecting the correct boundary
        % conditions for each question.
        function obj = ChemicalReaction(m,n)
            % If question 2, initialise values including the size, S.
            if m == 2
                obj.X = 0;
                obj.Y = 0;
                obj.S = ChemicalReaction.S0(n);
                obj.n_max = 10^8;
                obj.t_max = 10;
                obj.t_delta = 0.005;
                obj.c = [4800*obj.S, 20, 3e-5/obj.S^2, 6];
                % Else, question 1, intialise values including cdX, cdY, x0
                % and y0 for respective parts (a) and (b)
            else
                obj.X = ChemicalReaction.x0(n);
                obj.Y = ChemicalReaction.y0(n);
                obj.t_max = 4;
                obj.t_delta = 0.01;
                obj.n_max =  10^4;
                obj.cdx = ChemicalReaction.CDX(n);
                obj.cdy = ChemicalReaction.CDY(n);
            end
            % Initialising common values for all question parts.
            obj.t = ChemicalReaction.t_min;
            obj.count = 0;
            obj.t_sample = 0;
            obj.n_sample = 0;
            % Creating a random number stream, using Mersenne Twister random
            % number generator algorithm, and constant seed for
            % reporoducibility in producing different realisations.
            stream = RandStream('mt19937ar','Seed',ChemicalReaction.seed);
            RandStream.setGlobalStream(stream);
        end
        % Simulation section which calculates and stores X and Y population
        % values with their corresponding sample times.
        function obj = calculate(obj,m)
            % loop while count < n_max and t < t_max
            while (obj.count < obj.n_max && obj.t < obj.t_max)
                % If the property S is empty in the object, it must be for
                % question 1, so question 1 propensities are calculated.
                if isempty(obj.S) == 1
                    % h function in propensity is inputted for X and Y
                    obj.h = [0.5*obj.X*(obj.X-1) , obj.Y];
                    % multiplied by cdX and cdY to form the propensities
                    obj.a = [obj.cdx obj.cdy].*obj.h;
                    % If S has been assigned a value, object corresponds to
                    % question 2.
                else
                    % h function in propensity is inputted for X and Y
                    obj.h = [1,obj.X,0.5*obj.X*(obj.X-1)...
                        *obj.Y,obj.X];
                    % multiplied by cdX and cdY to form the propensities
                    obj.a = obj.c.*obj.h;
                end
                % total propensity is calculated
                obj.a0 = sum(obj.a);
                % two random numbers [0,1] are generated
                obj.random = rand(2,1);
                % the random time till next reaction is calculated
                obj.tau = log(1/obj.random(1))/obj.a0;
                %cumulative sum of a returns a 1 x 2 array of the
                %propensity of X and the propensity of X and Y summed
                %together. The inequality returns a logial 1 x 2 array
                %showing 0 if r2*a0 is bigger than each cumuative
                %propensity. If the sum of these logical values is 2 then
                %r2*a0 is smaller than both summations in
                %the inequality equation and the first reaction
                %occurs. If the sum is 1 then the inequality holds and the
                %second reaction occurs. Similarly for question 2 except mu
                %values can take {1,2,3,4} as there are 4 reactions. if mu
                %= 0, no reaction occurs.
                obj.mu = sum(obj.random(2)*obj.a0 <= cumsum(obj.a));
                % t is updated by tau
                obj.t = obj.t + obj.tau;
                % loop while t > t_sample and until the first t_sample
                % value over t_max is recorded as other values are
                % disregarded.
                while (obj.t > obj.t_sample && (obj.t_sample - obj.t_delta) < obj.t_max)
                    % n_sample is updated
                    obj.n_sample = obj.n_sample + 1;
                    % Y and X populations are stored alongside the sample
                    % time
                    obj.ydata(obj.n_sample,:) = obj.Y;
                    obj.xdata(obj.n_sample,:) = obj.X;
                    obj.tdata(obj.n_sample,:) = obj.t_sample;
                    % the mean of the recorded data with each iteration is
                    % calculated and stored in an array
                    obj.mean_(obj.n_sample,:) = sum(obj.xdata)/obj.n_sample;
                    % if question 1, also find the predicted values for X and Y for
                    % each sample time and save in an array for NRMSD
                    % calculation
                    if isempty(obj.S) == 1
                        obj.X_pred(obj.n_sample,:) = ((1/ChemicalReaction.x0(m))+obj.cdx.*obj.tdata(obj.n_sample)).^(-1);
                        obj.Y_pred(obj.n_sample,:) = ChemicalReaction.y0(m).*exp(-obj.cdy.*obj.tdata(obj.n_sample));
                        % if question 2, calculated the cumulative standard
                        % deviation and skewness values to be plotted in
                        % question 2(c).
                    else
                        obj.standarddeviation(obj.n_sample,:) = sqrt(sum((obj.xdata - obj.mean_(obj.n_sample,:)).^2)/(obj.n_sample-1));
                        obj.skewness(obj.n_sample,:) = (obj.n_sample*(sum(((obj.xdata-obj.mean_(obj.n_sample,:))/obj.standarddeviation(obj.n_sample,:)).^3)))/((obj.n_sample-1)*(obj.n_sample-2));
                    end
                    % update the sample time by t_delta.
                    obj.t_sample = obj.t_sample + obj.t_delta;
                end
                % if question 1, perform reaction depending on value of mu
                % if mu is 0, nothing occurs.
                if isempty(obj.S) == 1
                    switch obj.mu
                        case 2 % (reaction 1)
                            % X is reduced by 2
                            obj.X =obj.X - 2;
                        case 1 % (reaction 2)
                            % Y is reduced by 1
                            obj.Y = obj.Y - 1;
                    end
                    % if question 2, perform reaction depending on value of
                    % mu. If mu is 0, nothing occurs
                else
                    switch obj.mu
                        case 4  % (reaction 1)
                            % X is increased by 1
                            obj.X = obj.X + 1;
                        case 3 % (reaction 2)
                            % X is reduced by 1
                            obj.X = obj.X - 1;
                            % Y is increased by 1
                            obj.Y = obj.Y + 1;
                        case 2 % (reaction 3)
                            % X is increased by 1
                            obj.X = obj.X + 1;
                            % Y is reduced by 1
                            obj.Y = obj.Y - 1;
                        case 1 % (reaction 4)
                            % X is reduced by 1
                            obj.X = obj.X - 1;
                    end
                end
                % the number of time steps is updated
                obj.count = obj.count+1;
            end
        end
        % calculating NRMSD values and plotting the realisations in
        % question 1.
        function [nrmsd_y,nrmsd_x] = realisation1(obj,i)
            % calculating NRMSD for X and Y realisations
            nrmsd_y = (sqrt((sum((obj.Y_pred - obj.ydata).^2))/obj.n_sample))/obj.mean_(obj.n_sample);
            nrmsd_x = (sqrt((sum((obj.X_pred - obj.xdata).^2))/obj.n_sample))/obj.mean_(obj.n_sample);
            % plotting the Y molecule realisation with a stair-step plot as
            % values are discrete
            stairs(obj.tdata,obj.ydata,'color',ChemicalReaction.blues(i,:))
            ylabel('Number of X and Y molecules','Interpreter','latex');
            xlabel('Time Units','interpreter','latex')
            xlim([0 obj.t_max])
            hold on
            % plotting the X molecule realisation in the same plot with a stair-step plot as
            % values are discrete
            stairs(obj.tdata, obj.xdata,'color',ChemicalReaction.reds(i,:))
            hold on
        end
        % function to calculate and plot the determinsitic solutions
        function obj = deterministic(obj,m)
            % calculating deterministic plot as an array
            obj.X_det = ((1/ChemicalReaction.x0(m))+obj.cdx.*ChemicalReaction.t_det).^(-1);
            obj.Y_det = ChemicalReaction.y0(m).*exp(-obj.cdy.*ChemicalReaction.t_det);
            % plotting both X and Y deterministic solutions
            det = plot(ChemicalReaction.t_det,obj.Y_det, ChemicalReaction.t_det, obj.X_det,'linewidth',1.5);
            set(det, {'color'}, {[1 0 1]; [0 1 0]});
            legend({'Number of Y Molecules','Number of X molecules','Deterministic Solution for Y','Deterministic Solution for X'},'Interpreter','latex');
        end
        % plotting realisations in question 2
        function obj = realisation2(obj)
            % plot Y realisation in a stair-step plot
            stair_y = stairs(obj.tdata,obj.ydata,'color',(1/256).*[0 142 204]);
            ylabel('Number of X and Y molecules','Interpreter','latex')
            xlabel('Time Units','interpreter','latex')
            % dynamically title each plot depending on the size value
            title(strcat('S =  ',num2str(obj.S)),'interpreter','latex')
            xlim([0 obj.t_max])
            hold on
            % plot X realisation in a stair-step plot
            stair_x = stairs(obj.tdata,obj.xdata,'color',(1/256).*[210 31 60]);
            legend([stair_y stair_x],{'Y Molecules','X Molecules'},'Interpreter','latex','Location','southoutside')
        end
        % plotting the probability distribution from the data as well as
        % generated normal and kernel distribution and performing a
        % chi-squared test for normality
        function obj = prob_dist(obj)
            % plot data in a normalised histogram for actual probability
            % distribution
            p1 = histogram(obj.xdata,'normalization','probability','BinWidth',1,'EdgeColor','none','Facecolor',(1/256).*[0 142 204]);
            xlabel('Population of Molecule X','interpreter','latex')
            ylabel('Probability Density','interpreter','latex')
            hold on
            % fitting the data into a normal distribution and outputting
            % the mean, standard deviation and their respective confidence
            % intervals
            norm = fitdist(obj.xdata,'Normal')
            % creating a normal probability distribution function
            ynorm = pdf(norm,obj.xdata);
            % plotting the normal probability distribution
            p2 = line(obj.xdata,ynorm,'linestyle','-.','color','r');
            hold on
            % fitting the data into a kernel distribution
            ker = fitdist(obj.xdata,'Kernel');
            % creating a kernel probability distribution
            yker = pdf(ker,obj.xdata);
            % plotting the kernel distribution function
            p3 = line(obj.xdata,yker,'linestyle','--','color',[1 0 1]);
            legend([p1 p2 p3],{'X Population','Normal','Kernel'},'interpreter','latex')
            xlim([0 400])
            % chi-squared goodness of fit test for a normal distribution
            % with a 5 percent significance level
            [h0,p,stats] = chi2gof(obj.xdata,'CDF',norm,'Alpha',0.05,'NBins',400)
        end
        % plotting the estimated mean, estimated standard deviation and
        % skewness
        function obj = quantities(obj)
            % finding the length of the X population data
            [m,~] = size(obj.xdata);
            % plotting quantities in one figure as sub plots for
            % conciseness.
            subplot(1,3,1)
            plot(obj.tdata,obj.mean_,'color',(1/256).*[0 142 204]);
            ylabel('Cumulative Mean of X Population','interpreter','latex')
            xlabel('Time Units','interpreter','latex')
            xlim([0 obj.t_max])
            subplot(1,3,2)
            plot(obj.tdata,obj.standarddeviation,'color',(1/256).*[0 142 204]);
            ylabel('Cumulative Standard Deviation of X Population','interpreter','latex')
            xlabel('Time Units','interpreter','latex')
            xlim([0 obj.t_max])
            subplot(1,3,3)
            plot(obj.tdata,obj.skewness,'color',(1/256).*[0 142 204]);
            ylabel('Cumulative Skewness of X Population','interpreter','latex')
            xlabel('Time Units','interpreter','latex')
            xlim([0 obj.t_max])
            % outputting the estimated mean, estimated standard deviation
            % and skewness values for the whole dataset
            disp(['Mean = ',num2str(obj.mean_(m)), ' Standard Deviation = ',num2str(obj.standarddeviation(m)), ...
                ' Skewness = ', num2str(obj.skewness(m))]);
        end
        % Calculating the proportion of the data above and below the mean, within 1, 2 and 3 standard deviations of the mean and the median are calculated.
        function obj = additional_quantities(obj)
            a_mean = (sum(obj.xdata > obj.mean_(length(obj.xdata)))/length(obj.xdata))*100;
            b_mean = (sum(obj.xdata < obj.mean_(length(obj.xdata)))/length(obj.xdata))*100;
            std1 = 100 - (sum(obj.xdata < (obj.mean_(length(obj.xdata))-obj.standarddeviation(length(obj.xdata))))/length(obj.xdata))*100 - (sum(obj.xdata > (obj.mean_(length(obj.xdata))+obj.standarddeviation(length(obj.xdata))))/length(obj.xdata))*100;
            std2 = 100 - (sum(obj.xdata < (obj.mean_(length(obj.xdata))-2*(obj.standarddeviation(length(obj.xdata)))))/length(obj.xdata))*100 - (sum(obj.xdata > (obj.mean_(length(obj.xdata))+2*(obj.standarddeviation(length(obj.xdata)))))/length(obj.xdata))*100;
            std3 = 100 - (sum(obj.xdata < (obj.mean_(length(obj.xdata))-3*(obj.standarddeviation(length(obj.xdata)))))/length(obj.xdata))*100 - (sum(obj.xdata > (obj.mean_(length(obj.xdata))+3*(obj.standarddeviation(length(obj.xdata)))))/length(obj.xdata))*100;
            med = median(obj.xdata);
            disp([' Proportion of data above mean = ',num2str(a_mean),' Proportion of data below mean = ',num2str(b_mean), ' Proportion of data within 1 standard deviation = ',num2str(std1), ...
                 ' Proportion of data within 2 standard deviations = ',num2str(std2),' Proportion of data within 3 standard deviations = ',num2str(std3),'Median = ',num2str(med)]);
        end
        % calculating and plotting the power spectral density of the longer
        % simulation
        function obj = psd(obj,k)
            % subtracting the mean from the data so the mean is now zero so
            % natural spike at 0 frequency is accounted for
            obj.xdatanorm = obj.xdata - mean(obj.xdata);
            obj.ydatanorm = obj.ydata - mean(obj.ydata);
            % calculating the power spectral density estimates and
            % frequency vectors for X and Y
            [obj.pxx_x,obj.fx] = periodogram(obj.xdatanorm,[],[],obj.n_sample/obj.t_max);
            [obj.pxx_y,obj.fy] = periodogram(obj.ydatanorm,[],[],obj.n_sample/obj.t_max);
            % plotting the periodogram for all S values in one figure as
            % subplots.
            subplot(1,3,k);
            p1 = plot(obj.fy,obj.pxx_y,'color',(1/256).*[0 142 204]);
            hold on
            p2 = plot(obj.fx,obj.pxx_x,'color',(1/256).*[210 31 60]);
            title(strcat('S =  ',num2str(ChemicalReaction.S0(k))),'interpreter','latex')
            hold on
            xlim([0 6]);
            legend([p1 p2],{'Y','X'},'interpreter','latex','location','northeast')
            % adding one x-axis and y-axis title only
            if k == 1
                ylabel('PSD Magnitude (power / frequency)','interpreter','latex')
            elseif k == 2
                xlabel('Frequency (cycles / time unit)','interpreter','latex')
            end
        end
    end
    % calculating the mean NRMSD for X and Y in question 1(a) and (b) of
    % all realisations does not need to be a property in an object it is a
    % static method
    methods (Static)
        % calculating mean NMRSD for all X and Y realisations in questions
        % 1(a) and 1(b)
        function nrmsd = mean_nrmsd(ax,ay,bx,by)
            % storing as a 1 x 4 array of values
            nrmsd = [sum(ay)/length(ay),sum(ax)/length(ax), ...
                sum(by)/length(by),sum(bx)/length(bx)];
        end
    end
end