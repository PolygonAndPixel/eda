( find ./ -name '*.cpp' -print0 | xargs -0 cat ) | wc -l
( find ./ -name '*.h' -print0 | xargs -0 cat ) | wc -l
( find ./ -name '*.xml' -print0 | xargs -0 cat ) | wc -l