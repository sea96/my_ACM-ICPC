"Linux vim ~/.vimrc
set cin nu et ts=4 sw=4 sts=4 noswapfile nobackup cursorline
set background=dark
syntax on
map<F4> :exec"!g++ -std=c++11 -o2 % -o %<"<cr>
map<F5> :exec"!./%<"<cr>

"Windows gVim
set go= cin nu ts=4 sw=4 sts=4 et noswapfile nobackup cursorline
set bs=eol,start,indent
set fileencodings=utf-8,ucs-bom,cp936,gb18030
set guifont=Consolas:h12
syntax on
colorscheme molokai
map <F4> :exec "!g++ -std=c++11 -o2 % -o %<"<cr>
map <F5> :exec "!%<"<cr>
map <F2> :!javac %<cr>
map <F3> :!java %<<cr>
