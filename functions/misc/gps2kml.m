function gps2kml(filename,gpstime,Lat,Lon,Alt,style,opt11,opt12,opt21,opt22)
%
%example:
%[time,lat,lon,height]=textread('gpsOut.dat','%n%n%n%n');
%gps2kml('gpsOut',time,lat,lon,height,'o-r','MarkerSize',10,'LineWidth',3);

% function gps2kml(filename,Lat,Lon,style,opt11,opt12,opt21,opt22)
%
% Description: creates a file in kmz format that can be opened into Google Earth.
%    GEplot uses the same syntax as the traditional plot function but
%    requires Latitude and Longitudd (WGS84) instead of x and y (note that Lat is
%    the first argument).
%    If you need to convert from UTM to Lat/Lon you may use utm2deg.m also
%    available at Matlab Central
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Argument checking
%
error(nargchk(5, 10, nargin));  %5 arguments required, 10 maximum
n1=length(Lat);
n2=length(Lon);
n3=length(Alt);
if (n1~=n2)
   error('Lat and Lon vectors should have the same length');
end

if (n1~=n3)
    error ('Altitude should have the same length');
end

if (nargin==7 || nargin==9)
    error('size arguments must be "MarkerSize" or "LineWidth" strings followed by a number');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% symbol size and line width
%
markersize=7; %matlab default
linewidth=2;  %matlab default is 0.5, too thin for map overlay
if (nargin==8)
    if (strcmpi(opt11,'markersize')==1)
        markersize=opt12;
    elseif (strcmpi(opt11,'linewidth')==1)
        linewidth=opt12;
    else
        error('size arguments must be "MarkerSize" or "LineWidth" strings followed by a number');
    end
end
if (nargin==10)
    if (strcmpi(opt21,'markersize')==1)
        markersize=opt22;
    elseif (strcmpi(opt21,'linewidth')==1)
        linewidth=opt22;
    else
        error('size arguments must be "MarkerSize" or "LineWidth" strings followed by a number');
    end
end
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% symbol, line style and color
%
symbol='none';
iconfilename='none';
linestyle='-';
color='b';
colorstring='ffff0000'; 

if (nargin>=6)
    %linestyle
    if (strfind(style,'-'))
        linestyle='-';
    else
        linestyle='none';
    end

    %symbol
    if (strfind(style,'.')), symbol='.'; iconfilename='dot'; end
    if (strfind(style,'o')), symbol='o'; iconfilename='circle'; end
    if (strfind(style,'x')), symbol='x'; iconfilename='x'; end
    if (strfind(style,'+')), symbol='+'; iconfilename='plus'; end
    if (strfind(style,'*')), symbol='*'; iconfilename='star'; end
    if (strfind(style,'s')), symbol='s'; iconfilename='square'; end
    if (strfind(style,'d')), symbol='d'; iconfilename='diamond'; end
    if (strfind(style,'S')), symbol='S'; iconfilename='Ssquare'; end
    if (strfind(style,'D')), symbol='D'; iconfilename='Sdiamon'; end
    if (strfind(style,'O')), symbol='O'; iconfilename='dot'; end
    if (strfind(style,'0')), symbol='O'; iconfilename='dot'; end

    %color
    if (strfind(style,'b')), color='b'; colorstring='ffff0000'; end
    if (strfind(style,'g')), color='g'; colorstring='ff00ff00'; end
    if (strfind(style,'r')), color='r'; colorstring='ff0000ff'; end
    if (strfind(style,'c')), color='c'; colorstring='ffffff00'; end
    if (strfind(style,'m')), color='m'; colorstring='ffff00ff'; end
    if (strfind(style,'y')), color='y'; colorstring='ff00ffff'; end
    if (strfind(style,'k')), color='k'; colorstring='ff000000'; end
    if (strfind(style,'w')), color='w'; colorstring='ffffffff'; end
end



iconfilename='http://maps.google.com/mapfiles/kml/pal5/icon6.png';
if (symbol=='.') 
    markersize=markersize/5;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Creating kml file
%
fp=fopen(strcat(filename,'.kml'),'w');
if (fp==-1)
    message=sprint('Unable to open file %s.kml',filename);
    error(message);
end
fprintf(fp,'<?xml version="1.0" encoding="UTF-8"?>\n');

fprintf(fp,'<kml xmlns="http://www.opengis.net/kml/2.2">\n');
% fprintf(fp,'<kml xmlns="http://earth.google.com/kml/2.1">\n');

fprintf(fp,'<Document>\n');
fprintf(fp,'<name>%s</name>\n',strcat(filename,'.kml'));
fprintf(fp,'<description>Convert from GPS to KML (Marvell)</description>\n');
%
%Symbol styles definition
fprintf(fp,'<Style id="mystyle">\n');

fprintf(fp,'   <IconStyle>\n');
fprintf(fp,'    <color>%s</color>\n', colorstring);
%fprintf(fp,'     <colorMode>random</colorMode>\n');
fprintf(fp,'      <scale>%6.3f</scale>\n',markersize/20); %scale adjusted for .png image sizes
fprintf(fp,'      <Icon><href>%s</href></Icon>\n',iconfilename);      
fprintf(fp,'   </IconStyle>\n');


fprintf(fp,'   <LineStyle>\n');
fprintf(fp,'      <color>%s</color>\n',colorstring);
fprintf(fp,'      <width>%d</width>\n',linewidth);
fprintf(fp,'   </LineStyle>\n');


fprintf(fp,'    <LabelStyle>\n');
fprintf(fp,'    <color>%s</color>\n', colorstring);
fprintf(fp,'    <scale> 0.3 </scale>\n');
fprintf(fp,'    </LabelStyle>\n');

fprintf(fp,'</Style>\n');
fprintf(fp,'\n');


% if (linestyle=='-')
%     fprintf(fp,'    <Placemark>\n');
%     fprintf(fp,'      <description><![CDATA[ Field Test Trajectory]]>\n  We shall not cease from exploration \n </description>\n');
%     fprintf(fp,'      <name>Test Trajectory</name>\n');
%     fprintf(fp,'      <visibility>1</visibility>\n');
%     fprintf(fp,'      <open>1</open>\n');
%     fprintf(fp,'      <styleUrl>mystyle</styleUrl>\n');
%     fprintf(fp,'      <LineString>\n');
%     fprintf(fp,'        <extrude>1</extrude>\n');
%     fprintf(fp,'        <tessellate>0</tessellate>\n');
%     fprintf(fp,'        <altitudeMode>clampToGround</altitudeMode>\n');
%     fprintf(fp,'        <coordinates>\n');
%     for k=1:n1
%       fprintf(fp,'%15.8f, %15.8f, %15.8f\n',Lon(k),Lat(k),Alt(k));
%     end
%     fprintf(fp,'        </coordinates>\n');
%     fprintf(fp,'      </LineString>\n');
%     fprintf(fp,'    </Placemark>\n');
% end


if (strcmp(symbol,'none')==0)
fprintf(fp,'    <Folder>\n');
fprintf(fp,'      <name>Solution points</name>\n');

for k=1:n1
    fprintf(fp,'      <Placemark>\n');
    fprintf(fp,'         <description><![CDATA[Location Info\nGPSTime=%15.6f s\nLon = %15.6f deg\nLat = %15.6f deg\nAlt = %15.6f m \n\n \n]]></description>\n',gpstime(k),Lon(k),Lat(k),Alt(k));
    
%     Mode = GPS<br/>SOW = 174298.16<br/>RTC = 2765.59 s<br/>UTC = 2009-08-18 00:24:43.16<br/></description>
% 			<TimeStamp><when>2009-08-18T00:24:43.16Z</when></TimeStamp>
   

    fprintf(fp,'         <name>P %d</name>\n',k);  %you may add point labels here
    
    fprintf(fp,'<TimeStamp>\n <when>%12.3f\n</when></TimeStamp>\n',gpstime(k));
    
    fprintf(fp,'         <visibility>1</visibility>\n');
    fprintf(fp,'         <open>1</open>\n');
    fprintf(fp,'         <styleUrl>#mystyle</styleUrl>\n');
    fprintf(fp,'         <Point>\n');
    fprintf(fp,'           <coordinates>\n');
    fprintf(fp,'%15.8f, %15.8f, %15.8f\n',Lon(k),Lat(k),Alt(k));
    fprintf(fp,'           </coordinates>\n');
    fprintf(fp,'         </Point>\n');
    fprintf(fp,'      </Placemark>\n');
end
fprintf(fp,'    </Folder>\n');
end

fprintf(fp,'</Document>\n');
fprintf(fp,'</kml>\n');

fclose(fp);

if (strcmp(symbol,'none')==1)
   zip(filename,{strcat(filename,'.kml')});
% else
%    zip(filename,{strcat(filename,'.kml'), iconfilename});
% zip(filename,{strcat(filename,'.kml')});
end
% movefile(strcat(filename,'.zip'),strcat(filename,'.kmz'));
% delete(strcat(filename,'.kml'));
