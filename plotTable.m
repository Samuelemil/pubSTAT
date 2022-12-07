function plotTable(t)
% Function that plots tables as either HTML if publishing or pretty table
% format if not
%  Author: Samuel E Schmidt sschmidt@hst.aau.dk

%%
if ispublishing()
    disp('<html></html>')

    disp('<html>  </html>')
    disp([   strjoin(table2html(t)) ])
    disp('<html></html>')
    disp('<html> </html>')
else

    for i=1:size(t,2)

        t.(i)=categorical(t.(i));
        t.(i)(isundefined(t.(i)))={'-  -  -  -'};
    end

    disp(t)

end

end