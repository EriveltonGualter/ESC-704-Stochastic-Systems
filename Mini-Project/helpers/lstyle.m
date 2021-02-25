%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [line,marker]=lstyle(s)
% decode the linestyles

line=[];marker=[];
markers='+o*.xsd^v><ph';
if length(s)>3, return;end; % just in case a property was passed

for k=1:length(markers);
   if ~isempty(findstr(s,markers(k)))
      marker=markers(k);
      break
   end
end
s=[s,'   '];
if ~isempty(findstr(s,'--'))
   line='--';
elseif ~isempty(findstr(s,'-.'))
   line='-.';
   if (length(findstr(s,'.'))<2 & marker=='.'), marker=[];end;
elseif ~isempty(findstr(s,':'))
   line=':';
elseif ~isempty(findstr(s,'-'))
   line='-'; 
end