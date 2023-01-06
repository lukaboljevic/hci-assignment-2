% From EEGMMI DB:
%
% Tasks:
%     Task 1 (open and close left or right fist)
%     Task 2 (imagine opening and closing left or right fist)
%     Task 3 (open and close both fists or both feet)
%     Task 4 (imagine opening and closing both fists or both feet)
%
% T0, T1, T2 intervals:
%     T0 corresponds to rest
%     T1 corresponds to onset of motion (real or imagined) of
%         the left fist (in runs 3, 4, 7, 8, 11, and 12)
%         both fists (in runs 5, 6, 9, 10, 13, and 14)
%     T2 corresponds to onset of motion (real or imagined) of
%         the right fist (in runs 3, 4, 7, 8, 11, and 12)
%         both feet (in runs 5, 6, 9, 10, 13, and 14)


subject = "S053";
readWhichList = [3 4];  % 3 - task 1, 4 - task 2, 5 - task 3, 6 - task 4
nList = [20 35 50 60 75 90 100 115];
plot = false;

for i=1:length(readWhichList)
    readWhich = readWhichList(i);
    for j=1:length(nList)
        n = nList(j);
        [featVecFile, refFile] = features(subject, readWhich, "data", n, plot);
        classifyActivities(featVecFile, refFile, {1, 1}, 10, 30, 0, plot);
    end
end

% featVecFile = strcat("results/", subject, "/", subject, "-task", num2str(readWhich-2), "-", num2str(n), "-featv.txt");
% refFile = strcat("results/", subject, "/", subject, "-task", num2str(readWhich-2), "-", num2str(n), "-ref.txt");

