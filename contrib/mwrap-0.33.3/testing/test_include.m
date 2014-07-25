function test_include()

mex_id_ = 'o int = add2(i int)';
[j] = test_includemex(mex_id_, 2);
tassert(j == 4, 'Include test');

% ================================================================
function tassert(pred, msg)

if ~pred, fprintf('Failure: %s\n', msg); end
