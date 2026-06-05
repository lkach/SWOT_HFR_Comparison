function cascade_figures(N_f)

% Pos1 = get(gcf,'Position');

disp('This is configured for an extended screen at negative x coordinates.')

if length(N_f) == 1
    for ii = 1:N_f
        figure(ii)
        GCF_N = gcf;
        GCF_N.Position(1) = -1900; % 100 for one-screen
        GCF_N.Position(2) =  1100;
        set(GCF_N,'Position',GCF_N.Position + [(1500/N_f)*(ii-1) -(500/N_f)*(ii-1) 0 0])
    end
    
else
    jj = 1;
    for ii = N_f
        figure(ii)
        GCF_N = gcf;
        GCF_N.Position(1) = -1900; % 100 for one-screen
        GCF_N.Position(2) =  1100;
        set(GCF_N,'Position',GCF_N.Position + [(1500/length(N_f))*(jj-1) -(500/length(N_f))*(jj-1) 0 0])
        jj = jj + 1;
    end
end

end