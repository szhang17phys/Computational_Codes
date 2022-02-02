set et
set sw=2
set smarttab
set smartindent
set bg=dark
:nmap <C-N><C-N> :set invnumber <CR>
set nohlsearch
syntax on
map <C-k> <C-w><Up>
map <C-j> <C-w><Down>
map <C-l> <C-w><Right>
map <C-h> <C-w><Left>
if has("autocmd")
  au BufReadPost * if line("'\"") > 0 && line("'\"") <= line("$")
    \| exe "normal! g'\"" | endif
endif
autocmd FileType * setlocal formatoptions-=c formatoptions-=r formatoptions-=o
