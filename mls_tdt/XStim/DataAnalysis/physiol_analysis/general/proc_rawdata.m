clear
tests = [361 364 367 370 381 382];
rec_filename = '899r02.rec';

for test_num = 1:length(tests)
   proc_raw(rec_filename, tests(test_num));
end