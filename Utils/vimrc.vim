'vim ~/.vimrc
set cin nu et ts=4 sw=4 sts=4 noswapfile nobackup cursorline
set background=dark
syntax on
map<F4> :exec"!g++ -std=c++11 -o2 % -o %<"<cr>
map<F5> :exec"!./%<"<cr>
