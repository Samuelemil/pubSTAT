function out = ispublishing()
  st  = dbstack;
  out = any(arrayfun(@(x)strcmp(x.file,'publish.p'),st)); 
end