

%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Multi-Objetective Atomic Orbital Search (MOAOS)            %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2022a (MacOs)               %
%                                                                   %
%                      Author and programmer                        %
%                ---------------------------------                  %
%                      Nima Khodadadi (ʘ‿ʘ)                         %
%                             e-Mail                                %
%                ---------------------------------                  %
%                         nkhod002@fiu.edu                          %
%                                                                   %
%                            Homepage                               %
%                ---------------------------------                  %
%                    https://nimakhodadadi.com                      %
%                                                                   %
%                                                                   %
%                                                                   %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% ----------------------------------------------------------------------- %





clear;
close all;
clc;
% format long g

% Initial parameters of the MOAOS algorithm
for j=8
    max_iter=100;
    Pop=100;
    ArchiveMaxSize=100;

    nvar=30;
    %     obj_no=2;

    method=3;
    if method==3

        Archive_F1=load(sprintf('P%d.txt',j));
    end
    obj_no=2;
    mkdir (sprintf('P%d',j))
    %-------------------------- MOAOS-----------------------------------------
    for i=1% Numbver of independent runs

        [Archive_F]=MAOS(max_iter,ArchiveMaxSize,Pop,nvar,method,j);
        if numel(Archive_F')==2
            continue
        end
        Archive_F=Archive_F';


        if obj_no==2
            plot(Archive_F1(:,1),Archive_F1(:,2),'Color','g','LineWidth',4);
            hold on
            plot(Archive_F(:,1),Archive_F(:,2),'ro','LineWidth',1,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','r',...
                'Marker','o',...
                'MarkerSize',10);
            legend('True PF','Obtained PF');
            title(sprintf('MOAOS FOR MMF%d PROBLEM',j));
            xlabel('obj_1');
            ylabel('obj_2');

            hold off
        end
        if obj_no==3
            plot3(Archive_F1(:,1),Archive_F1(:,2),Archive_F1(:,3),'Color','g','LineWidth',1);
            hold on
            plot3(Archive_F(:,1),Archive_F(:,2),Archive_F(:,3),'ro','LineWidth',1,...
                'MarkerEdgeColor','b',...
                'MarkerFaceColor','r',...
                'Marker','o',...
                'MarkerSize',10);
            legend('True PF','Obtained PF');
            title(sprintf('MOAOS FOR MMF%d PROBLEM',j));
            xlabel('obj_1');
            ylabel('obj_2');
            zlabel('obj_3');

            hold off
        end


        savefig(sprintf('P%d/fig_%d.fig',j,i));

    end




    save(sprintf('P%d/result_P%d.mat',j,j));

end