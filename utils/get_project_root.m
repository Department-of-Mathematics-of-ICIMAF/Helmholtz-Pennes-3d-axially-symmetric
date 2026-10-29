function root_dir = get_project_root()
%GET_PROJECT_ROOT Devuelve la raiz del repositorio.
    root_dir = fileparts(fileparts(mfilename('fullpath')));
end
