function save_figure(fig, output_dir, base_name)
%SAVE_FIGURE Guarda una figura en formatos FIG y PNG, y despues la cierra.
    if ~ishandle(fig)
        error('save_figure:InvalidFigure', 'La figura indicada no es valida.');
    end
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    saveas(fig, fullfile(output_dir, [base_name, '.fig']));
    print(fig, fullfile(output_dir, [base_name, '.png']), '-dpng', '-r150');
    close(fig);
end
