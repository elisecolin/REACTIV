%  (c) Elise Koeniguer 
%     « Copyright (c) 2016, Elise Koeniguer (Onera) » 
%     This file is part of REACTIV.
% 
%     REACTIV is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     REACTIV is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with GeFolki in the file copying.txt.  
%     If not, see <http://www.gnu.org/licenses/gpl.txt>.
%
% ------------------------------------------------------------------------
%
%     REACTIV est un logiciel libre ; vous pouvez le redistribuer ou le
%     modifier suivant les termes de la GNU General Public License telle
%     que publiée par la Free Software Foundation ; soit la version 3 de la
%     licence, soit (à votre gré) toute version ultérieure.
% 
%     REACTIV est distribué dans l'espoir qu'il sera utile, mais SANS
%     AUCUNE GARANTIE ; sans même la garantie tacite de QUALITÉ MARCHANDE
%     ou d'ADÉQUATION à UN BUT PARTICULIER. Consultez la GNU General Public
%     License pour plus de détails.
% 
%     Vous devez avoir reçu une copie de la GNU General Public License en
%     même temps que Gefolki dans le fichier copying.txt ; si ce n'est pas
%     le cas, consultez % <http://www.gnu.org/licenses/gpl.txt

% Define the directory for your images
Rep='/c7/MEDUSE/IDF/Saclay/Sentinel/Pile_From_SNAP/SLC/PileRecalee/DESCENDING/IW1/Recalage-Elise/AmplitudesTif/';
% Select List of Files using command "dir"
List=dir([Rep,'*VH_amplitude.tif']);


% Select the dates 
N=length(List);
dateN=zeros(N,1);
for k=1:N
    y=str2num(List(k).name(1:4));
    m=str2num(List(k).name(5:6));
    d=str2num(List(k).name(7:8));
    dateN(k)=datenum(y,m,d);
end

D=zeros(N,6);
D(:,1)=y;
D(:,2)=m;
D(:,3)=d;

C=reactiv_list(Rep,List,0,0,dateN);



