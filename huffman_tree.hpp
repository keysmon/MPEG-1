#include <iostream>
#include <fstream>
#include <stack>
#include <queue>
#include <string.h>
#include <cstdio>
//#include <windows.h>
using namespace std;

struct node{
    int num;
    unsigned char uch;
    node *left;
    node *right;
    node *parent;
    node() {};

    node(int n,unsigned char uc, node *l, node *r, node *p) : num(n), uch(uc), left(l), right(r), parent(p) {}
};

class huf_tree{
public:
    node *root;
    string st;
    
	huf_tree() {
	    node * n = new node(0, 255, NULL, NULL, NULL);
	    root = n;
	    st = "";
	}
	
	bool swap(node* &a,node* &b) {
	    if(a == b || a == NULL || b == NULL || b->parent == a || a->parent == b) {
	    	return false;
		}
	    int pos1, pos2;
	    if(b->parent->left == b) {
	        pos1 = 1;
	    }
	    else pos1 = 0;
	    if(a->parent->left == a) {
	        pos2 = 1;
	    }
	    else pos2 = 0;
	    node *p = new node;
	    p = a->parent;
	    a->parent = b->parent;
	    b->parent = p;
	    if(pos1 == 1) {
	        a->parent->left = a;
	    }
	    else a->parent->right = a;
	    if(pos2 == 1) {
	        b->parent->left = b;
	    }
	    else b->parent->right = b;

	    return true;
	}
	
	string addnode(unsigned char c) { 
	    if(root->left == NULL) {
	        node *l = new node(0, 250, NULL, NULL, root);
	        root->left = l;
	        node *r = new node(1, c, NULL, NULL, root);
	        root->right = r;
	        this->st += c;
	        string s = "        ";
	        for(int i = 7; i >= 0; i --) {
	            s[i] = '0' + (c & 1);
	            c /= 2;
	        }
	        return s;
	    }
	    if(is_exist(c)) {
	        node *n = getnode(c, 1);
	        string r = CreateCode(n);
	        rebuild_tree(n);
	        return r;
	    }
	    else{//ÐÂ×Ö·û 
	        node *n = getnode(250, 0);
	        string s = "        ";
	        unsigned char c1 = c;
	        for(int i = 7; i >= 0; i--) {
	            s[i] = (c1 % 2) + '0';
	            c1 /= 2;
	        }
	        s = CreateCode(n) + s;
	        node *l = new node(0, 250, NULL, NULL, n);
	        n->left = l;
	        node *r = new node(1, c, NULL, NULL, n);
	        n->right = r;
	        this->st += c;
	        rebuild_tree(n);
	        return s;
	    }
	}
	
	void rebuild_tree(node *n) {
	    while(n != root) {
	        stack<node*> s;
	        queue<node*> q;
	        q.push(this->root);
	        while(!q.empty()) {
                node* tmp = q.front();
                s.push(tmp);
                if(tmp->right != NULL) {
                    q.push(tmp->right);
                }
                if(tmp->left != NULL) {
                    q.push(tmp->left);
                }
				q.pop();
	        }
	        node* tmp1 = s.top();
	        s.pop();
	        while(!s.empty()) {
	            tmp1 = s.top();
	            if(tmp1 == n) break;
	            s.pop();
	        }
	        node* nswap = s.top();
	        while(s.size()>1) {
	            if(tmp1->num != s.top()->num) break;
	            nswap = s.top();
	            s.pop();
	        }
	        swap(nswap, tmp1);
	        tmp1->num++;
	        n = tmp1->parent;
	    }
	}
	
	
	node *getnode(unsigned char c, int flag) {
	    queue<node*> q;
	    q.push(root);
	    while(!q.empty()) {
			node* t = q.front();
			if(flag == 1) { 
				if(t->uch == c) {
	            	return t;
	            }
			}
			else {
				if(t->uch == c && (t->left == NULL) && (t->right == NULL)) {
	                return t;
	            }
			}
            
            if(t->right!=NULL) q.push(t->right);
            if(t->left!=NULL) q.push(t->left);
            q.pop();
	    }
	    return NULL;
	}
	
	
	string CreateCode(node* p) {
	    string s = "", path = "";
	    for ( p; p != this->root; p = p->parent) {
	        if(p->parent->left == p) {
	            s += '0';
	        }
	        else s += '1';
	    }
	    for(int i = s.size()-1; i >= 0; i--) {
	        path += s[i];
	    }
	    return path;
	}
	
	
	bool is_exist(unsigned char c) {
	    for(int i = 0; i < st.size(); i ++) {
	        if(c == st[i]) return true;
	    }
	    return false;
	}
};

unsigned char chTobit(char s[8]) {
    char result = 0;
    for(int i = 0; i < 8; i++) {
        result = result *2;
        result += s[i] - 48;
    }
    return result;
}


