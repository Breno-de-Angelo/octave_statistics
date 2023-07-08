
myfile = fopen ("dados.txt", "r");

i=0;
wspd = [];
wvht = [];
while ~feof(myfile)
    i++;
    line = fgetl(myfile);
    if(i <= 2)
      continue;
    endif
    linhaFormatada = textscan(line, "%f");
    valoresLinha = cell2mat(linhaFormatada);
    valorWSPD = valoresLinha(7);
    valorWVHT = valoresLinha(9);
    wspd = [wspd; valorWSPD];
    wvht = [wvht; valorWVHT];
end

save("-binary","vetorWSPD", "wspd");
save("-binary", "vetorWVHT", "wvht");

fclose(myfile);
