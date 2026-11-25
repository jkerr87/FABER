function filename = export_FABER_vehicle_report(out,templateFile,filename)
% Refresh one orbit tab in a copy of the MA_DV_Vehicle Excel layout.
arguments
    out (1,1) struct
    templateFile (1,1) string
    filename (1,1) string
end
assert(isfile(templateFile),'Provide an existing vehicle report template.');
[parent,name,extension] = fileparts(filename);
assert(strlength(name)>0 && strcmpi(extension,'.xlsx'), ...
    'Provide a separate .xlsx output path.');
if isempty(parent) || strlength(string(parent))==0, parent = pwd; end
if ~isfolder(parent), mkdir(parent); end
[~,folderInfo] = fileattrib(parent);
filename = string(fullfile(folderInfo.Name,name+extension));
[~,templateInfo] = fileattrib(templateFile);
assert(~strcmpi(filename,templateInfo.Name),'Output must differ from the original template.');

% Keep the other orbit tab if it was generated on a previous run. Stage
% edits beside the destination so a failed refresh leaves the prior copy.
staged = string(tempname(parent))+".xlsx";
cleanup = onCleanup(@() removeStagedFile(staged)); %#ok<NASGU>
existing = isfile(filename);
source = templateFile;
if existing, source = filename; end
[ok,message] = copyfile(source,staged);
assert(ok,'%s',message);

row = [4:17 29:35 52:54];
label = ["Injection","Client Changes","Departure RPO from Atlante", ...
    "Far field rendezvous to client","RPO to client","Far field rendezvous to Atlante", ...
    "RPO to Atlante","Station keeping","Relative station keeping","Safety","Abort","CAM", ...
    "Momentum Unloading","Disposal","Injection","Client Changes","Station keeping", ...
    "Relative station keeping","CAM","Momentum Unloading","Disposal", ...
    "Assembly relocations to host structure","Assembly relocations to Atlante","Momenutm Unloading"];
vehicle = [repmat("Taxi",1,14),repmat("ATLANTE",1,7),repmat("Kiki",1,3)];
maneuver = ["Initial orbit correction","Client orbit transfer","Atlante Departure", ...
    "Far-field rendezvous to Client","RPO to Client","Far-field rendezvous to ATLANTE", ...
    "RPO to ATLANTE","Absolute stationkeeping","Relative stationkeeping", ...
    "Safety contingency and recovery","Aborted approach","Collision avoidance", ...
    "Momentum unloading","Disposal","Initial orbit correction","Client orbit transfer", ...
    "Absolute stationkeeping","Relative stationkeeping","Collision avoidance", ...
    "Momentum unloading","Disposal","Loaded relocation","Unloaded return","Momentum unloading"];

sheets = ["ISS-Like","SSO-Like"];
whichSheet = 1+(out.cfg.i_target_deg>90);
for j = 1:numel(sheets)
    if existing && j~=whichSheet, continue; end
    actual = readcell(staged,'Sheet',sheets(j),'Range','A4:A54');
    assert(all(strcmpi(strtrim(string(actual(row-3))),label(:))), ...
        'Unexpected vehicle report layout in the %s sheet.',sheets(j));
end

B = print_FABER_budget(out,"");
values = cell(numel(row),5);
for j = 1:numel(row)
    match = B.Vehicle==vehicle(j) & B.Maneuver==maneuver(j);
    assert(nnz(match)==1,'Expected one budget row for %s / %s.',vehicle(j),maneuver(j));
    b = B(match,:);
    time = b.ModeledTime_min;
    if maneuver(j)=="Disposal" && isfield(out,'disposal')
        selected = out.disposal.Selected & out.disposal.Vehicle==vehicle(j);
        if nnz(selected)==1, time = sprintf('%.4f yr',out.disposal.Decay_years(selected)); end
    end
    values(j,:) = {b.RawDV_mps,b.MarginedDV_mps-b.RawDV_mps,b.MarginedDV_mps,time,b.Notes};
end

for j = 1:numel(sheets)
    if existing && j~=whichSheet, continue; end
    writeReportCells({upper(sheets(j))},staged,sheets(j),'D1');
    content = values;
    if j~=whichSheet
        % A new copy starts with an empty results area on the unrun tab.
        content = repmat({NaN},numel(row),5);
    end
    for range = {1:14,15:21,22:24}
        indices = range{1};
        writeReportCells(content(indices,:),staged,sheets(j), ...
            sprintf('E%d:I%d',row(indices(1)),row(indices(end))));
    end
    writeReportCells({"End-of-life disposal to the modeled decay cutoff"}, ...
        staged,sheets(j),'B17');
end
if isfield(out,'transferFiring')
    writetable(out.transferFiring,staged,'Sheet',"Transfer Firing "+extractBefore(sheets(whichSheet),'-'), ...
        'WriteMode','overwritesheet');
end
if isfield(out,'rpoChecks') && isfield(out.rpoChecks,'propulsion')
    writetable(out.rpoChecks.propulsion,staged,'Sheet',"RPO Propulsion "+extractBefore(sheets(whichSheet),'-'), ...
        'WriteMode','overwritesheet');
end
[ok,message] = movefile(staged,filename,'f');
assert(ok,'%s',message);
fprintf('Vehicle report updated: %s (%s)\n',filename,sheets(whichSheet));
end

function writeReportCells(values,filename,sheet,range)
writecell(values,filename,'Sheet',sheet,'Range',range, ...
    'PreserveFormat',true,'AutoFitWidth',false);
end

function removeStagedFile(path)
if isfile(path), delete(path); end
end
