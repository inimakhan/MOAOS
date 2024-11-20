
%__________________________________________________________________ %
%                          Multi-Objective                          %
%        Multi-Objetective Atomic Orbital Search (MOAOS)            %
%                                                                   %
%                                                                   %
%                  Developed in MATLAB R2021a (MacOs)               %
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



% ----------------------------------------------------------------------- %
function [Index SubIndex]=GetGridIndex(particle,G)

    c=particle.Cost;
    
    nobj=numel(c);
    ngrid=numel(G(1).Upper);
    
    str=['sub2ind(' mat2str(ones(1,nobj)*ngrid)];

    SubIndex=zeros(1,nobj);
    for j=1:nobj
        
        U=G(j).Upper;
        
        i=find(c(j)<U,1,'first');
        
        SubIndex(j)=i;
        
        str=[str ',' num2str(i)];
    end
    
    str=[str ');'];
    
    Index=eval(str);
    
end