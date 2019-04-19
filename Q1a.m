%Author: Brennan Powers

close all
clear 

disp(' ')
disp('Enjoy localization with my particle filter')
disp(' ')

% number of particles
particleNumber = 1000;

fid = fopen('localization.log', 'r');
data = textscan(fid,'%s %f %f %f');
fclose(fid);

landmark_spotted = 0;
rng(0,'twister'); % makes randoms consistent for testing
min_range = -2500;
max_range = 2500;
r(:,1) = (max_range - min_range).*rand(particleNumber,1) + min_range;
r(:,2) = (max_range - min_range).*rand(particleNumber,1) + min_range;
min_radians = 0;
max_radians = 2*pi;
r(:,3) = (max_radians - min_radians).*rand(particleNumber,1) + min_radians;
% so col 1 of r is x coord, col 2 is y coord, col 3 is orientation in rad.

hold on;
%axis([-2500 2500 -2500 2500]);
axis([-3500 3500 -3500 3500]);
current_plot = scatter(r(:,1),r(:,2), '.', 'b');
plot(-1500,2250,'g*'); % landmark 0
plot(1500,2250,'g*');  % landmark 1
plot(1500,0,'g*');     % landmark 2
plot(-1500,0,'g*');    % landmark 3
%axis('manual');
%axis([-2500 2500 -2500 2500]);
%hold off;

landmark0 = [-1500 2250];
landmark1 = [1500 2250];
landmark2 = [1500 0];
landmark3 = [-1500 0];

for i=1 : size(data{1,1},1)
    pause(0.05);
    if data{1,1}(i,1) == "turn"
        min_radians = -.05;
        max_radians = .05;
        r(:,3) = mod( (((max_radians - min_radians).*rand(particleNumber,1) + min_radians) ... 
            + (r(:,3) + data{1,2}(i,1))),(2*pi));
        
        %r(:,3) = mod((r(:,3) + data{1,2}(i,1)),(2*pi));
    elseif data{1,1}(i,1) == "step"
        min_range = -10;
        max_range = 10;
        r(:,1) = (((max_range - min_range).*rand(particleNumber,1) + min_range)) ...
            + (r(:,1) + ( data{1,2}(i,1) * cos(r(:,3))) );
        r(:,2) = (((max_range - min_range).*rand(particleNumber,1) + min_range)) ... 
            + (r(:,2) + ( data{1,2}(i,1) * sin(r(:,3))) );
        
        %r(:,1) = r(:,1) + ( data{1,2}(i,1) * cos(r(:,3)) );
        %r(:,2) = r(:,2) + ( data{1,2}(i,1) * sin(r(:,3)) );
    elseif data{1,1}(i,1) == "landmark"
        landmark_spotted = 1;
        if data{1,2}(i,1) == 0

            r(:,4) = pdist2( [ r(:,1) r(:,2) ], [landmark0(1,1),landmark0(1,2)]);
            
            
            r(:,5) = atan((r(:,2)-landmark0(1,2))./(r(:,1)-landmark0(1,1))) ...
                - (pi * (r(:,1) > landmark0(1,1)) );
            %r(:,6) = rad2deg(r(:,5));
            
            %r(:,7) = ((r(:,4) >= data{1,3}(i,1)*0.85) & (r(:,4) <= data{1,3}(i,1)*1.15)) ... 
            %    + ((r(:,3)+data{1,4}(i,1)+.05 <= data{1,4}(i,1)) & (r(:,3)+data{1,4}(i,1)-.05 >= data{1,4}(i,1)));
                        
%             x = -1500;
%             y = 2250;
%             radius = data{1,3}(i,1);
%             %---------------
%             % plots circles
%             th = 0:pi/50:2*pi;
%             xunit = radius * cos(th) + x;
%             yunit = radius * sin(th) + y;
%             h = plot(xunit, yunit, 'r');
%             pause(1);
%             %---------------
        elseif data{1,2}(i,1) == 1

            r(:,4) = pdist2( [ r(:,1) r(:,2) ], [landmark1(1,1),landmark1(1,2)]);
            r(:,5) = atan((r(:,2)-landmark1(1,2))./(r(:,1)-landmark1(1,1))) ...
                - (pi * (r(:,1) > landmark1(1,1)) );
           
            %r(:,6) = rad2deg(r(:,5));
            
            %r(:,7) = ((r(:,4) >= data{1,3}(i,1)*0.85) & (r(:,4) <= data{1,3}(i,1)*1.15)) ... 
            %    + ((r(:,3)+data{1,4}(i,1)+.05 <= data{1,4}(i,1)) & (r(:,3)+data{1,4}(i,1)-.05 >= data{1,4}(i,1)));

%             x = 1500;
%             y = 2250;
%             radius = data{1,3}(i,1);
%             %---------------
%             % plots circles
%             th = 0:pi/50:2*pi;
%             xunit = radius * cos(th) + x;
%             yunit = radius * sin(th) + y;
%             h = plot(xunit, yunit);
%             %---------------
%             pause(1);
        elseif data{1,2}(i,1) == 2

            r(:,4) = pdist2( [ r(:,1) r(:,2) ], [landmark2(1,1),landmark2(1,2)]);
            r(:,5) = atan((r(:,2)-landmark2(1,2))./(r(:,1)-landmark2(1,1))) ...
                - (pi * (r(:,1) > landmark2(1,1)) );
            
            %r(:,6) = rad2deg(r(:,5));
            
            %r(:,7) = ((r(:,4) >= data{1,3}(i,1)*0.85) & (r(:,4) <= data{1,3}(i,1)*1.15)) ... 
            %    + ((r(:,3)+data{1,4}(i,1)+.05 <= data{1,4}(i,1)) & (r(:,3)+data{1,4}(i,1)-.05 >= data{1,4}(i,1)));
            
%             x = 1500;
%             y = 0;
%             radius = data{1,3}(i,1);
%             %---------------
%             % plots circles
%             th = 0:pi/50:2*pi;
%             xunit = radius * cos(th) + x;
%             yunit = radius * sin(th) + y;
%             h = plot(xunit, yunit);
%             %---------------
%             pause(1);
        elseif data{1,2}(i,1) == 3

            r(:,4) = pdist2( [ r(:,1) r(:,2) ], [landmark3(1,1),landmark3(1,2)]);
            r(:,5) = atan((r(:,2)-landmark3(1,2))./(r(:,1)-landmark3(1,1))) ...
                - (pi * (r(:,1) > landmark3(1,1)) );
            
            %r(:,6) = rad2deg(r(:,5));
            
            %r(:,7) = ((r(:,4) >= data{1,3}(i,1)*0.85) & (r(:,4) <= data{1,3}(i,1)*1.15)) ... 
            %    + ((r(:,3)+data{1,4}(i,1)+.05 <= data{1,4}(i,1)) & (r(:,3)+data{1,4}(i,1)-.05 >= data{1,4}(i,1)));

%             x = -1500;
%             y = 0;
%             radius = data{1,3}(i,1);
%             %---------------
%             % plots circles
%             th = 0:pi/50:2*pi;
%             xunit = radius * cos(th) + x;
%             yunit = radius * sin(th) + y;
%             h = plot(xunit, yunit);
%             %---------------
%             pause(1);
        else
            disp('Some error in reading the file.');
        end
        
    else
        %return;
    end
    
    if landmark_spotted == 0
        %set(current_plot, 'Visible', 'off');
        hold off;
        %axis('manual');
        %axis([-2500 2500 -2500 2500]);
        %hold on;
        current_plot = scatter(r(:,1),r(:,2), '.', 'blue');
        axis([-3500 3500 -3500 3500]);
        %axis([-2500 2500 -2500 2500]);
        hold on;
        %hold on;
        plot(-1500,2250,'g*');
        plot(1500,2250,'g*');
        plot(1500,0,'g*');
        plot(-1500,0,'g*');
        %hold on;
    else
%         hold off;
%         current_plot = scatter(r(:,1),r(:,2), '.', 'blue');
%         axis([-2500 2500 -2500 2500]);
%         hold on;
%         plot(-1500,2250,'g*');
%         plot(1500,2250,'g*');
%         plot(1500,0,'g*');
%         plot(-1500,0,'g*');
%         scatter(r(:,1).*r(:,7),r(:,2).*r(:,7), 'o', 'red');
        

        z = [(data{1,3}(i,1)-r(:,4)) , (data{1,4}(i,1)-r(:,5))];
        p = mvnpdf(z, [0 0], [1.15*data{1,3}(i,1) 0; 0 .05]);
        %y(:,1) = randsample(r(:,1),1000,true,p);
        %y(:,2) = randsample(r(:,2),1000,true,p);
        %y(:,3) = randsample(r(:,3),1000,true,p);
        if max(p) == 0
            p(:,1) = .00001;
        end
        
        min_range = -10;
        max_range = 10;
        r(:,1) = ((max_range - min_range).*rand(particleNumber,1) + min_range) ...
            + randsample(r(:,1),particleNumber,true,p);
        r(:,2) = ((max_range - min_range).*rand(particleNumber,1) + min_range) ... 
            + randsample(r(:,2),particleNumber,true,p);
        min_radians = -.05;
        max_radians = .05;
        r(:,3) = mod( (((max_radians - min_radians).*rand(particleNumber,1) + min_radians) ... 
            + randsample(r(:,3),particleNumber,true,p)),(2*pi));

        hold off;
        current_plot = scatter(r(:,1),r(:,2), '.', 'blue');
        %axis([-2500 2500 -2500 2500]);
        axis([-3500 3500 -3500 3500]);
        hold on;
        plot(-1500,2250,'g*');
        plot(1500,2250,'g*');
        plot(1500,0,'g*');
        plot(-1500,0,'g*');
        
        landmark_spotted = 0;
        
    end
    
    
    
end



